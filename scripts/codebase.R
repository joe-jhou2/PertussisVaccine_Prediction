
###############################################################################
#' Load packages and required directory paths
###############################################################################

#install.packages("pacman")
library(pacman)

p_load(devtools, tidyverse, Hmisc, BiocManager)
p_load(corrplot, ggpubr, impute, Biobase, limma, mice)
p_load_gh("mani2012/BatchQC")

## Batch correction
p_load(vsn, lme4, sva)

## MOfa
p_load(MOFA2, psych)

## modelling
p_load(omicade4, mogsa, RSpectra, lubridate, glmnet)


### Source code from MCIA model
source('https://raw.githubusercontent.com/akonstodata/mcia_mbpca/main/R/MCIA_mbpca_extra.R')

## Read all master paths
dir_rds_objects <- paste0(base_dir, "objects/")

## Read all housekeeping functions


###############################################################################
#' A function for correcting date format
###############################################################################
parse_mdy_2digit_years <- function(x) {
  library(stringr)
  library(lubridate)
  
  x <- str_trim(x)
  dates <- suppressWarnings(mdy(x))
  # Fix 2068+ to 1968-style
  dates <- if_else(!is.na(dates) & year(dates) > Sys.Date() %>% year(),
                   dates %m-% years(100),
                   dates)
  return(dates)
}

###############################################################################
# Define a function for processing Subject data
###############################################################################
process_subject_data <- function(subject_data) {
  subject_data %>%
    mutate(
      age_at_boost = year(date_of_boost) - year(year_of_birth),
      age = age_at_boost - (month(date_of_boost) < month(year_of_birth) | 
                              (month(date_of_boost) == month(year_of_birth) & day(date_of_boost) < day(year_of_birth)))
    )
}

###############################################################################
#' A function for principal variance component analysis
#'
#' The function is written based on the 'pvcaBatchAssess' function of the PVCA R package
#' and slightly changed to make it more efficient and flexible for sequencing read counts data.
#' (http://watson.nci.nih.gov/bioc_mirror/packages/release/bioc/manuals/pvca/man/pvca.pdf)
#'
#' @param counts The Normalized(e.g. TMM)/ log-transformed reads count matrix from sequencing data (row:gene/feature, col:sample)
#' @param meta  The Meta data matrix containing predictor variables (row:sample, col:predictor)
#' @param threshold The proportion of the variation in read counts explained by top k PCs. This value determines the number of top PCs to be used in pvca.
#' @param inter TRUE/FALSE - include/do not include pairwise interactions of predictors
#'
#' @return std.prop.val The vector of proportions of variation explained by each predictor.
#'
#' @export
#'
###############################################################################

pvca_run <- function (abatch, batch.factors, threshold) 
{
  library(vsn)
  library(lme4)
  
  theDataMatrixCentered <- t(apply(abatch, 1, scale, center=TRUE, scale=FALSE))
  
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues/eigenValuesSum
  expInfo <- pData(abatch)[, batch.factors]
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1) {
    my_sum_2 = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold) {
      my_counter_2 = my_counter_2 + 1
    }
  }
  if (my_counter_2 < 3) {
    pc_n = 3
  }
  else {
    pc_n = my_counter_2
  }
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
                                               pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                                                        i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                  ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                   expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                  na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, 
                                                           j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, 
                                                                 j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames))
}

pvca_run_rnaseq <- function (abatch, batch.factors, threshold) 
{
  library(vsn)
  library(lme4)
  
  theDataMatrix <- exprs(vsn2(abatch, verbose = FALSE, minDataPointsPerStratum = 15L))
  dataRowN <- nrow(theDataMatrix)
  dataColN <- ncol(theDataMatrix)
  theDataMatrixCentered <- matrix(data = 0, nrow = dataRowN, 
                                  ncol = dataColN)
  theDataMatrixCentered_transposed = apply(theDataMatrix, 1, scale, center = TRUE, scale = FALSE)
  theDataMatrixCentered = t(theDataMatrixCentered_transposed)
  
  #theDataMatrixCentered <- t(apply(abatch, 1, scale, center=TRUE, scale=FALSE))
  
  theDataCor <- cor(theDataMatrixCentered)
  eigenData <- eigen(theDataCor)
  eigenValues = eigenData$values
  ev_n <- length(eigenValues)
  eigenVectorsMatrix = eigenData$vectors
  eigenValuesSum = sum(eigenValues)
  percents_PCs = eigenValues/eigenValuesSum
  expInfo <- pData(abatch)[, batch.factors]
  exp_design <- as.data.frame(expInfo)
  expDesignRowN <- nrow(exp_design)
  expDesignColN <- ncol(exp_design)
  my_counter_2 = 0
  my_sum_2 = 1
  for (i in ev_n:1) {
    my_sum_2 = my_sum_2 - percents_PCs[i]
    if ((my_sum_2) <= threshold) {
      my_counter_2 = my_counter_2 + 1
    }
  }
  if (my_counter_2 < 3) {
    pc_n = 3
  }
  else {
    pc_n = my_counter_2
  }
  pc_data_matrix <- matrix(data = 0, nrow = (expDesignRowN * 
                                               pc_n), ncol = 1)
  mycounter = 0
  for (i in 1:pc_n) {
    for (j in 1:expDesignRowN) {
      mycounter <- mycounter + 1
      pc_data_matrix[mycounter, 1] = eigenVectorsMatrix[j, 
                                                        i]
    }
  }
  AAA <- exp_design[rep(1:expDesignRowN, pc_n), ]
  Data <- cbind(AAA, pc_data_matrix)
  variables <- c(colnames(exp_design))
  for (i in 1:length(variables)) {
    Data$variables[i] <- as.factor(Data$variables[i])
  }
  op <- options(warn = (-1))
  effects_n = expDesignColN + choose(expDesignColN, 2) + 1
  
  randomEffectsMatrix <- matrix(data = 0, nrow = pc_n, ncol = effects_n)
  model.func <- c()
  index <- 1
  for (i in 1:length(variables)) {
    mod = paste("(1|", variables[i], ")", sep = "")
    model.func[index] = mod
    index = index + 1
  }
  for (i in 1:(length(variables) - 1)) {
    for (j in (i + 1):length(variables)) {
      mod = paste("(1|", variables[i], ":", variables[j], 
                  ")", sep = "")
      model.func[index] = mod
      index = index + 1
    }
  }
  function.mods <- paste(model.func, collapse = " + ")
  for (i in 1:pc_n) {
    y = (((i - 1) * expDesignRowN) + 1)
    funct <- paste("pc_data_matrix", function.mods, sep = " ~ ")
    Rm1ML <- lmer(funct, Data[y:(((i - 1) * expDesignRowN) + 
                                   expDesignRowN), ], REML = TRUE, verbose = FALSE, 
                  na.action = na.omit)
    randomEffects <- Rm1ML
    randomEffectsMatrix[i, ] <- c(unlist(VarCorr(Rm1ML)), 
                                  resid = sigma(Rm1ML)^2)
  }
  effectsNames <- c(names(getME(Rm1ML, "cnms")), "resid")
  randomEffectsMatrixStdze <- matrix(data = 0, nrow = pc_n, 
                                     ncol = effects_n)
  for (i in 1:pc_n) {
    mySum = sum(randomEffectsMatrix[i, ])
    for (j in 1:effects_n) {
      randomEffectsMatrixStdze[i, j] = randomEffectsMatrix[i, 
                                                           j]/mySum
    }
  }
  randomEffectsMatrixWtProp <- matrix(data = 0, nrow = pc_n, 
                                      ncol = effects_n)
  for (i in 1:pc_n) {
    weight = eigenValues[i]/eigenValuesSum
    for (j in 1:effects_n) {
      randomEffectsMatrixWtProp[i, j] = randomEffectsMatrixStdze[i, 
                                                                 j] * weight
    }
  }
  randomEffectsSums <- matrix(data = 0, nrow = 1, ncol = effects_n)
  randomEffectsSums <- colSums(randomEffectsMatrixWtProp)
  totalSum = sum(randomEffectsSums)
  randomEffectsMatrixWtAveProp <- matrix(data = 0, nrow = 1, 
                                         ncol = effects_n)
  for (j in 1:effects_n) {
    randomEffectsMatrixWtAveProp[j] = randomEffectsSums[j]/totalSum
  }
  return(list(dat = randomEffectsMatrixWtAveProp, label = effectsNames))
}

plot_pvca <- function(pvcaObj, plot_title) {
  
  pvca_barplot_data <- data.frame(explained = pvcaObj$dat[1,],
                                  effect    = pvcaObj$label) %>%
    arrange(explained) %>%
    mutate(effect = factor(effect, levels = effect))
  
  pvca_barplot <- ggplot(data    = pvca_barplot_data,
                         mapping = aes(x = effect, y = explained)) +
    geom_bar(stat = "identity", fill="blue") +
    geom_text(aes(label = signif(explained, digits = 3)),
              nudge_y   = 0.01,
              size      = 4) +
    labs(x = "Effects", y = "Proportion of the variance explained") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1,
                                     size  = 12)) +
    ggtitle(plot_title)
  
  return(pvca_barplot)
  
}


plot_pca = function(count_data, metadata, plot_title = "PRCOMP Analysis"){
  
  pca_res <- prcomp(t(count_data), scale = TRUE)
  #fviz_eig(pca_res)
  
  var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
  
  pca_plot_PC12 <- pca_res$x %>% 
    as.data.frame %>%
    rownames_to_column("specimen_id") %>%
    mutate(specimen_id = as.numeric(specimen_id)) %>%
    left_join(metadata) %>%
    #separate(specimen_id,c("dataset")) %>%
    ggplot(aes(x=PC1,y=PC2)) + geom_point(aes(color=dataset),size=4) +
    theme_bw(base_size=16) + 
    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
    theme(legend.position="top") +
    ggtitle(plot_title)
  
  pca_plot_PC13 <- pca_res$x %>% 
    as.data.frame %>%
    rownames_to_column("specimen_id") %>%
    mutate(specimen_id = as.numeric(specimen_id)) %>%
    left_join(metadata) %>%
    #separate(specimen_id,c("dataset")) %>%
    ggplot(aes(x=PC1,y=PC3)) + geom_point(aes(color=dataset),size=4) +
    theme_bw(base_size=16) + 
    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
         y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
    theme(legend.position="top") +
    ggtitle(plot_title)
  
  pca_plot_PC23 <- pca_res$x %>% 
    as.data.frame %>%
    rownames_to_column("specimen_id") %>%
    mutate(specimen_id = as.numeric(specimen_id)) %>%
    left_join(metadata) %>%
    #separate(specimen_id,c("dataset")) %>%
    ggplot(aes(x=PC2,y=PC3)) + geom_point(aes(color=dataset),size=4) +
    theme_bw(base_size=16) + 
    labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
         y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
    theme(legend.position="top") +
    ggtitle(plot_title)
  
  return_plots = list(
    
    PC12 = pca_plot_PC12,
    PC13 = pca_plot_PC13,
    PC23 = pca_plot_PC23
    
  )
  
  return(return_plots)
  
}

pvca_analysis = function(mat_data, subject_specimen, batch.factors, plot_title = "PVCA analysis"){
  
  counts_data <- mat_data[rowMeans(is.na(mat_data)) < 1, ] %>%
    as.matrix() %>%
    impute.knn() %>%
    .$data
  
  subject_specimen <- subject_specimen %>%
    dplyr::filter(specimen_id %in% colnames(counts_data)) %>%
    dplyr::select(specimen_id, timepoint, infancy_vac, biological_sex, dataset) %>%
    dplyr::mutate(specimen_id1 = specimen_id) %>%
    column_to_rownames("specimen_id1") 
  
  subject_specimen <- subject_specimen[colnames(counts_data), ]
  
  phenoData = AnnotatedDataFrame(data = subject_specimen)
  
  exprset = ExpressionSet(as.matrix(counts_data), 
                          phenoData = phenoData)
  
  pvcaObj_after <- pvca_run(exprset, batch.factors = batch.factors, threshold =0.6) 
  
  print(plot_pvca(pvcaObj_after, paste0(plot_title)))
  
  
  #calculate principle components cell components
  pca_plots = plot_pca(counts_data, subject_specimen, paste0(plot_title))
  plot(pca_plots$PC12)
  plot(pca_plots$PC13)
  plot(pca_plots$PC23)
  
}

###############################################################################
#' Helper functions for Antibody levels data processing
#' to support PVCA analysis and batch correction etc
###############################################################################

processAbtiter = function(master_obj, BatchCorrection = TRUE){
  
  ## Assemble data
  count_data_long = master_obj$plasma_antibody_levels$long
  df_subject_specimen = master_obj$subject_specimen
  
  ## Reshape dataframe in wide format
  antibody_wide_normalized_pre <- count_data_long  %>%
    filter(specimen_id %in% df_subject_specimen$specimen_id) %>%
    dplyr::select(isotype_antigen, specimen_id, MFI_normalised) %>%
    pivot_wider(names_from = "isotype_antigen", values_from = MFI_normalised) %>%
    column_to_rownames("specimen_id")%>%
    t() 
  
  normalized_imputed = antibody_wide_normalized_pre[rowMeans(is.na(antibody_wide_normalized_pre)) < 1, ] %>%
    as.matrix() %>%
    impute.knn() %>%
    .$data
  
  if(BatchCorrection == TRUE){
    
    batch_lebels = as.data.frame(colnames(antibody_wide_normalized_pre)) %>%
      rename(specimen_id = starts_with("colnames")) %>%
      mutate(specimen_id = as.integer(specimen_id)) %>%
      left_join(df_subject_specimen) %>%
      dplyr::select(dataset)
    
    batchCorrected = ComBat(antibody_wide_normalized_pre, batch = batch_lebels$dataset)
    
    batchCorrected_imputed = batchCorrected[rowMeans(is.na(batchCorrected)) < 1, ] %>%
      as.matrix() %>%
      impute.knn() %>%
      .$data
    
    normalised_data = list(
      
      metadata = df_subject_specimen,
      raw_data = master_obj$plasma_antibody_levels$wide,
      normalized_data = normalized_imputed,
      batchCorrected_data = batchCorrected_imputed
    )
    
    return(normalised_data)
  }
  
  if(BatchCorrection == FALSE){
    
    ## Prepare return dataframe
    normalised_data = list(
      
      metadata = df_subject_specimen,
      raw_data = master_obj$plasma_antibody_levels$wide,
      normalized_data = antibody_normalized_imputed,
    )
    
    return(normalised_data)
  }
}

###############################################################################
#' Helper functions for Cell frequency data processing
#' to support data normalization, PVCA analysis and batch correction etc
###############################################################################

processCellFreq = function(master_obj, BatchCorrection = TRUE){
  
  ## Assemble data
  count_data_long = master_obj$pbmc_cell_frequency$long
  df_subject_specimen = master_obj$subject_specimen
  
  ## Perform median normalization
  cytof_median_D0 <- count_data_long %>%
    left_join(df_subject_specimen[c("specimen_id", "dataset")]) %>%
    filter(specimen_id %in% unique(df_subject_specimen[df_subject_specimen$planned_day_relative_to_boost == 0,]$specimen_id)) %>%
    group_by(dataset, cell_type_name)  %>%
    summarise(median = median(percent_live_cell, na.rm = T))
  
  cell_long_normalized_pre <-  count_data_long  %>%
    left_join(df_subject_specimen[c("specimen_id", "dataset")]) %>%
    left_join(cytof_median_D0) %>%
    mutate(percent_live_cell_normalized = if_else(is.na(percent_live_cell) == T, NA, percent_live_cell/median))
  
  ## Reshape dataframe in wide format
  cell_wide_normalized_pre <- cell_long_normalized_pre  %>%
    dplyr::select(cell_type_name, specimen_id, percent_live_cell_normalized) %>%
    pivot_wider(names_from = "cell_type_name", values_from = percent_live_cell_normalized) %>%
    column_to_rownames("specimen_id")%>%
    t() 
  
  cellFreq_normalized_imputed = cell_wide_normalized_pre[rowMeans(is.na(cell_wide_normalized_pre)) < 1, ] %>%
    as.matrix() %>%
    impute.knn() %>%
    .$data
  
  if(BatchCorrection == TRUE){
    
    batch_lebels = as.data.frame(colnames(cell_wide_normalized_pre)) %>%
      rename(specimen_id = starts_with("colnames")) %>%
      mutate(specimen_id = as.integer(specimen_id)) %>%
      left_join(df_subject_specimen) %>%
      dplyr::select(dataset)
    
    cellFreq_batchCorrected = ComBat(cell_wide_normalized_pre, batch = batch_lebels$dataset)
    
    cellFreq_batchCorrected_imputed = cellFreq_batchCorrected[rowMeans(is.na(cellFreq_batchCorrected)) < 1, ] %>%
      as.matrix() %>%
      impute.knn() %>%
      .$data
    
    normalised_data = list(
      
      metadata = df_subject_specimen,
      raw_data = master_obj$pbmc_cell_frequency$wide,
      normalized_data = cellFreq_normalized_imputed,
      batchCorrected_data = cellFreq_batchCorrected_imputed
    )
    
    return(normalised_data)
  }
  
  if(BatchCorrection == FALSE){
    
    ## Prepare return dataframe
    normalised_data = list(
      
      metadata = df_subject_specimen,
      raw_data = master_obj$pbmc_cell_frequency$wide,
      batchCorrected_data = cellFreq_normalized_imputed,
    )
    
    return(normalised_data)
  }
}

###############################################################################
#' Helper functions for Olink data processing
#' to support data normalization, PVCA analysis and batch correction etc
###############################################################################
processOlink <- function(master_obj, BatchCorrection = TRUE){
  
  ## Assemble data
  count_data_long = master_obj$plasma_cytokine_concentrations_by_olink$long
  df_subject_specimen = master_obj$subject_specimen
  
  ## Perform median normalization
  cytokine_median_D0 <- count_data_long %>%
    left_join(df_subject_specimen[c("specimen_id", "dataset")]) %>%
    filter(specimen_id %in% unique(df_subject_specimen[df_subject_specimen$planned_day_relative_to_boost == 0,]$specimen_id)) %>%
    group_by(dataset, protein_id)  %>%
    summarise(median = median(protein_expression, na.rm = T))
  
  cytokine_long_normalized_pre <-  count_data_long  %>%
    left_join(df_subject_specimen[c("specimen_id", "dataset")]) %>%
    left_join(cytokine_median_D0) %>%
    mutate(protein_expression_normalized = if_else(is.na(protein_expression) == T, NA, protein_expression/median))
  
  ## Reshape dataframe in wide format
  cytokine_wide_normalized_pre <- cytokine_long_normalized_pre  %>%
    dplyr::select(protein_id, specimen_id, protein_expression_normalized) %>%
    pivot_wider(names_from = "protein_id", values_from = protein_expression_normalized) %>%
    column_to_rownames("specimen_id")%>%
    t() 
  
  cytokineFreq_normalized_imputed = cytokine_wide_normalized_pre[rowMeans(is.na(cytokine_wide_normalized_pre)) < 1, ] %>%
    as.matrix() %>%
    impute.knn() %>%
    .$data
  
  if(BatchCorrection == TRUE){
    
    batch_lebels = as.data.frame(colnames(cytokine_wide_normalized_pre)) %>%
      rename(specimen_id = starts_with("colnames")) %>%
      mutate(specimen_id = as.integer(specimen_id)) %>%
      left_join(df_subject_specimen) %>%
      dplyr::select(dataset)
    
    cytokineFreq_batchCorrected = ComBat(cytokine_wide_normalized_pre, batch = batch_lebels$dataset)
    
    cytokineFreq_batchCorrected_imputed = cytokineFreq_batchCorrected[rowMeans(is.na(cytokineFreq_batchCorrected)) < 1, ] %>%
      as.matrix() %>%
      impute.knn() %>%
      .$data
    
    normalised_data = list(
      
      metadata = df_subject_specimen,
      raw_data = master_obj$pbmc_cytokine_frequency$wide,
      normalized_data = cytokineFreq_normalized_imputed,
      batchCorrected_data = cytokineFreq_batchCorrected_imputed
    )
    
    return(normalised_data)
  }
  
  if(BatchCorrection == FALSE){
    
    ## Prepare return dataframe
    normalised_data = list(
      
      metadata = df_subject_specimen,
      raw_data = master_obj$pbmc_cytokine_frequency$wide,
      normalized_data = cytokineFreq_normalized_imputed,
    )
    
    return(normalised_data)
  }
}

###############################################################################
#' Helper functions for Antibody titer data processing
#' to support PVCA analysis and batch correction etc
###############################################################################

pvca_analysis_rnaseq = function(mat_data, subject_specimen, batch.factors, plot_title = "PVCA analysis"){

  counts_data <- mat_data[rowMeans(is.na(mat_data)) < 1, ] %>%
    as.matrix() %>%
    impute.knn() %>%
    .$data
  
  subject_specimen <- subject_specimen %>%
    dplyr::filter(specimen_id %in% colnames(counts_data)) %>%
    dplyr::select(specimen_id, timepoint, infancy_vac, biological_sex, dataset) %>%
    dplyr::mutate(specimen_id1 = specimen_id) %>%
    column_to_rownames("specimen_id1") 
  
  subject_specimen <- subject_specimen[colnames(counts_data), ]
  
  phenoData = AnnotatedDataFrame(data = subject_specimen)
  
  exprset = ExpressionSet(as.matrix(counts_data), 
                          phenoData = phenoData)
  
  pvcaObj_after <- pvca_run_rnaseq(exprset, batch.factors = batch.factors, threshold =0.6) 
  
  print(plot_pvca(pvcaObj_after, paste0(plot_title)))
  
  #calculate principle components cell components
  pca_plots = plot_pca(counts_data, subject_specimen, paste0(plot_title))
  plot(pca_plots$PC12)
  plot(pca_plots$PC13)
  plot(pca_plots$PC23)
  
}

new_gs_pramod = function(data_input,mcia_out){
  
  nblocks=length(data_input)
  num_comps=dim(mcia_out$t)[2]
  blocks_out<-lapply(1:nblocks,matrix,data=NA,nrow=dim(data_input[[1]])[2],ncol=num_comps)
  names(blocks_out)<-names(data_input)
  for (i in 1:nblocks){
    name = names(data_input)[i]
    block = t(as.matrix(data.frame(data_input[name]))) %*% as.matrix(data.frame(mcia_out$pb[name]))
    blocks_out[[i]]<-block
    rownames(blocks_out[[i]])<-colnames(data_input[name])
  }
  weights = mcia_out$w
  
  gs = matrix(0,ncol(data_input$seq),num_comps)
  
  for (i in 1:num_comps){
    out_sum=0
    for (j in 1:nblocks){
      out = weights[j,i]*as.matrix(data.frame(blocks_out[[j]][,i]))
      out_sum=out_sum+out
    }
    gs[,i] = out_sum
  }
  return(gs)
}

subset_dataset = function(data_obj, datasets){
  
  metadata = data_obj$subject_specimen %>%
    filter(dataset %in% datasets)
  
  cytof_data_wide = data_obj$pbmc_cell_frequency$wide %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  cytof_data_long = data_obj$pbmc_cell_frequency$long %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  abtiter_data_wide = data_obj$plasma_antibody_levels$wide %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  abtiter_data_long = data_obj$plasma_antibody_levels$long %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  olink_data_wide = data_obj$plasma_cytokine_concentrations_by_olink$wide %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  olink_data_long = data_obj$plasma_cytokine_concentrations_by_olink$long %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  legendplex_data_wide = data_obj$plasma_cytokine_concentrations_by_legendplex$wide %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  legendplex_data_long = data_obj$plasma_cytokine_concentrations_by_legendplex$long %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  rnaseq_data_wide = data_obj$pbmc_gene_expression$wide %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  rnaseq_data_long = data_obj$pbmc_gene_expression$long %>%
    filter(specimen_id %in% metadata$specimen_id)
  
  subset_data <- list(
    
    subject_specimen = metadata,
    plasma_antibody_levels = list(
      
      wide = abtiter_data_wide,
      long = abtiter_data_long
    ),
    plasma_cytokine_concentrations_by_olink = list(
      
      wide = olink_data_wide,
      long = olink_data_long
    ),
    plasma_cytokine_concentrations_by_legendplex = list(
      
      wide = legendplex_data_wide,
      long = legendplex_data_long
    ),
    pbmc_cell_frequency_wide = list(
      
      wide = cytof_data_wide,
      long = cytof_data_long
    ),
    pbmc_gene_expression_wide = list(
      
      wide = rnaseq_data_wide,
      long = rnaseq_data_long
    )
  )
  
  return(subset_data)
  
}

mad_calculations = function(tpm_data, metadata, my_dataset, plot = TRUE){
  
  metadata_v2 <-  metadata %>% 
    filter(dataset %in% my_dataset,
           specimen_id %in% colnames(tpm_data)
    ) %>% 
    mutate(specimen_id = as.character(specimen_id))
  
  count_data <- tpm_data  %>% 
    dplyr::select(metadata_v2$specimen_id)
  
  #### Select top P fraction of genes with highest MAD
  P = 0.75
  mad_vals = apply(count_data, 1, mad)
  mad_cutoff = min(tail(sort(mad_vals), P*length(mad_vals)))
  df_mad <- data.frame(
    gene_id = rownames(count_data), 
    mean_exp = rowMeans((count_data)),
    mad_val = mad_vals)
  
  if(plot == TRUE){
    
    assay = paste0("RNAseq: ", my_dataset[1])
    p1 = ggplot(df_mad, aes(y=log1p(mad_val))) +
      geom_histogram(bins = 100) +
      coord_flip() + ggtitle(assay) + 
      geom_hline(yintercept = log1p(mad_cutoff), color="red")
    
    p2 = ggplot(df_mad, aes(log1p(mean_exp), log1p(mad_val))) + 
      #geom_point(size=0.5) + 
      geom_hex(bins=100) +
      geom_smooth() + ggtitle(assay) +
      geom_hline(yintercept = log1p(mad_cutoff), color="red") + 
      scale_fill_viridis_c(option = "A")
    
    print(ggarrange(p1, p2))
  }
  
  df_mad_return = df_mad %>%
    filter(mad_val >= mad_cutoff)
  
  return(df_mad_return)
}

# Recursive function to save data frames and explore lists
save_dataframes_to_tsv <- function(obj, prefix = "") {
  if (is.data.frame(obj)) {
    # If the object is a data frame, save it to a .tsv file
    write.table(obj, 
                file = paste0(dir_rds_objects, prefix, ".tsv"), 
                sep = "\t", 
                row.names = FALSE, 
                quote = FALSE)
  }
  else if (is.matrix(obj)) {
    # If the object is a data frame, save it to a .tsv file
    write.table(obj, 
                file = paste0(dir_rds_objects, prefix, ".tsv"), 
                sep = "\t", 
                row.names = TRUE, 
                quote = FALSE)
  }
  else if (is.list(obj)) {
    # If the object is a list, iterate over its elements
    for (name in names(obj)) {
      
      # Construct a new prefix for the file name based on the current name
      new_prefix <- ifelse(prefix == "", name, paste0(prefix, "_", name))
      
      print(paste0(prefix, "--->", new_prefix))
      # Recursive call with the element of the list
      save_dataframes_to_tsv(obj[[name]], new_prefix)
    }
  }
}

###############################################################################
# FUNCTIONS FOR MCIA_mbpca.Rmd
###############################################################################
library(corpcor)
library(RSpectra)

# Modified from: https://github.com/mengchen18/mogsa/blob/master/R/svd.nipals.R
# Description: Different methods are automatically selected to calculate svd
# according to matrix properties, such as missing values, object class, size. 

svd.solver <- function(
    x, nf, opts.svd.nipals=list(), opts.svds=list(), opts.fast.svd = list(), opts.svd=list()
) {
  if (miss.nf <- missing(nf))
    nf <- min(dim(x))
  
  smallnf <- nf < 3
  
  rt <- nrow(x)/ncol(x)
  if (any(is.na(x))) {
    r <- do.call(svd.nipals, c_mat(list(x=x, nf = nf), opts.svd.nipals))
    solver <- "svd.nipals"
    # cat("svd.nipals used. \n" )
  } else if (inherits(x, "Matrix") || smallnf) {
    r <- do.call(svds, c_mat(list(A=x, k= nf), opts.svds))
    solver <- "svds"
    # cat("svds used. \n" )
  } else if (length(x) >= 1e6 & abs(log10(rt)) >= 1) {
    if (!miss.nf)
      message("function 'fast.svd' used, all singular vectors will be computed")
    r <- do.call(fast.svd, c_mat(list(m=x), opts.fast.svd))
    solver <- "fast.svd"
    # cat("fast.svd used. \n" )
  } else {
    if (!miss.nf)
      message("Function 'svd' used, all singular vectors will be computed")
    r <- do.call(svd, c_mat(list(x=x), opts.svd))
    solver <- "svd"
    # cat("svd used. \n" )
  }
  if (inherits(x, "Matrix")) {
    r$u <- Matrix(r$u)
    r$v <- Matrix(r$v)
  }  
  attr(r, "solver") <- solver
  r
}

c_mat <- function(x, ...) {
  if (inherits(x, "Matrix"))
    as.numeric(x) else
      base::c(x, ...)
}

# modified from https://github.com/mengchen18/mogsa/blob/master/R/processOpt.R
# Description: data pre-processing for MCIA (within-block)
processOpt <-
  function(x, center=TRUE, scale=FALSE, num_comps, option = c("lambda1", "inertia", "uniform","lambda_all")) {
    
    opt <- match.arg(option)  
    
    if (is.null(names(x)))
      names(x) <- paste("data", 1:length(x), sep = "_")
    
    x <- lapply(x, function(x) {
      r <- scale(x, center = center, scale = scale)
      if (inherits(x, "Matrix"))
        r <- Matrix(r)
      r
    })
    if (opt == "lambda1") {
      w <- sapply(x, function(xx) 1/svd.solver(xx, nf = 1)$d)
    } else if (opt == "inertia") {
      w <- sapply(x, function(xx) 1/sqrt(sum(xx^2)))
    } else if (opt == "uniform") {
      w <- rep(1, length(x))
    } else if (opt == "lambda_all"){
      w <- sapply(x, function(xx) 1/sum(svd.solver(xx, nf = num_comps)$d))
    }
    mapply(SIMPLIFY = FALSE, function(xx, ww) xx*ww, xx=x, ww=w)
  }


# This is copied form made4 package in bioconductor
# v1.61

"array2ade4" <-
  function(dataset, pos=FALSE,  trans=FALSE){
    
    # Allows matrix, data.frame, ExpressionSet, marrayRaw to be read as data.frame
    if (!is.data.frame(dataset)) dataset<-getdata(dataset)
    
    if (any(is.na(dataset)))
      stop("Arraydata must not contain NA values. Use impute.knn in library(impute), KNNimpute from Troyanskaya et al., 2001 or LSimpute from Bo et al., 2004 to impute missing values\n")
    
    
    # COA needs table of positive data, will add real no to make +ve
    if(pos){
      if (any(dataset < 0)) {
        num<-round(min(dataset)-1)
        dataset<-dataset+abs(num)
      }
    }
    
    if(trans) {
      # Transpose matrix  (as BGA, CIA expects the samples to be in the rows)
      # dudi.nsc should not be transposed, use t.dudi instead to ensure row weight are equal
      # There is a horrible bug is dudi.pca/coa etc, if a dataset with vars>>cases is given
      # It can end abruptly crashing the session. This is a bug in sweep
      # There will now use t.dudi rather than transpose the data
      
      # using t convert data.frame to matrix and messes up affymetrix probe ID names
      # It changes all of the "-" to "." in probeids like AFFX-CreX-5_at
      # So save names and change the names of the transposed matrix
      
      colnam= colnames(dataset)
      rownam = rownames(dataset)               
      dataset<-t(dataset)		
      dimnames(dataset) = list(colnam, rownam) 
      
      # convert matrix to data.frame for ade4
      dataset <- as.data.frame(dataset)
      
      if (!is.data.frame(dataset)) stop("Problems checking dataset")
    }
    
    data.out<-dataset        
    return(data.out)
  }

# Custom functions

# Center and scale all feature x sample arrays in list 
center_scale<-function(x){
  data_out<-lapply(x, function(x) {
    r<-t(scale(t(x), center=TRUE, scale=TRUE))
  })
  return(data_out)}

# Preprocess arrays using nsc
nsc_prep<-function(data_input,num_comps){
  df.list<-data_input
  df.list <- lapply(df.list, as.data.frame)
  df.list <- lapply(df.list, array2ade4, pos = TRUE)
  coa.list <- lapply(df.list, dudi.nsc, scannf = FALSE, nf = num_comps)
  coa.list.t <- lapply(coa.list, ade4:::t.dudi)
  dfl <- lapply(coa.list, function(x) x$tab)
  return(dfl)
}

# Run MCIA or CPCA
# data_input is list of feature x samples arrays; samples must match
mcia_mbpca<-function(data_input,num_comps,preprocess='nsc',block_prep='inertia',
                     deflat_method='blockLoading'){
  if (preprocess=='nsc'){
    table_out<-nsc_prep(data_input, num_comps)
  }
  else{
    table_out<-center_scale(data_input)
  }
  # block-level preprocess 
  final_out<-processOpt(table_out,scale=FALSE,center=FALSE,num_comps,option=block_prep)
  
  mcia_out<-mbpca(final_out,ncomp=num_comps,k="all",method=deflat_method,
                  option="uniform",center=FALSE,scale=FALSE,
                  unit.p=TRUE,unit.obs=FALSE,moa=FALSE)
  rownames(mcia_out$t)<-colnames(data_input[[1]])
  mcia_out<-list(data_prep=final_out,mcia_result=mcia_out)
  return(mcia_out)
}

## Generate global scores from block loadings and block weights on new data
# data_input should have same features/omics as original
# data_input should be preprocessed in same way as in mcia_mbpca (akin to 'final_out')
new_gs<-function(data_input,mcia_out){
  nblocks=length(data_input)
  num_comps=dim(mcia_out$t)[2]
  #data_prep<-lapply(data_input,function(x) as.matrix(data.frame(x)))
  blocks_out<-lapply(1:nblocks,matrix,data=NA,nrow=dim(data_input[[1]])[2],ncol=num_comps)
  names(blocks_out)<-names(data_input)
  for (i in 1:nblocks){
    name = names(data_input)[i]
    block = t(as.matrix(data.frame(data_input[name]))) %*% as.matrix(data.frame(mcia_out$pb[name]))
    blocks_out[[i]]<-block
    rownames(blocks_out[[i]])<-colnames(data_input[name])
  }
  weights = mcia_out$w
  
  gs = matrix(0,nrow(mcia_out$t),ncol(mcia_out$t))
  for (i in 1:num_comps){
    out_sum=0
    for (j in 1:nblocks){
      out = weights[j,i]*as.matrix(data.frame(blocks_out[[j]][,i]))
      out_sum=out_sum+out
    }
    gs[,i] = out_sum
  }
  return(gs)
}