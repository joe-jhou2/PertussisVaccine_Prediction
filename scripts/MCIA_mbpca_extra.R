# FUNCTIONS FOR MCIA_mbpca.Rmd
# Anna Konstorum


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

