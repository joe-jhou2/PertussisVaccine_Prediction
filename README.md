# Pertussis Vaccine Response Prediction

## Overview

This project aims to build predictive models for immune response to Tdap vaccination using multi-modal datasets from acellular-Pertussis (aP) vs whole-cell Pertussis (wP) infancy-primed subjects. The analysis focuses on understanding vaccine efficacy and identifying potential biomarkers through machine learning approaches.

## Dataset Description

### Study Population
- **aP Group**: Individuals born after 1996 (acellular-Pertussis primed)
- **wP Group**: Individuals born before 1995 (whole-cell Pertussis primed)

### Study Design
- Baseline blood samples collected
- Tdap booster vaccine administered
- Follow-up blood samples at 1, 3, 7, 14, and 28 days post-vaccination

### Available Assays and Modalities

| Assay Type | Description | Method |
|------------|-------------|---------|
| Plasma Antibody Titers | Measured against Tdap at all time points | Luminex assay |
| Plasma Cytokine Concentrations | Cytokine analysis (Method 1) | Olink assay |
| Plasma Cytokine Concentrations | Cytokine analysis (Method 2) | Legendplex assay |
| PBMC Gene Expression | Bulk peripheral blood mononuclear cells | RNAseq analysis |
| PBMC Cell Frequency | Analysis of PBMC subset frequencies | Flow cytometry |
| T Cell Activation | T cell functional assessment | FluoroSpot assay |
| T Cell Polarization | T cell phenotyping | AIM assay |

## Prediction Targets

The project focuses on building predictive models for three key targets:

1. **Ab IgG_PT**: Plasma antibody titers against Tdap
2. **CCL3**: Gene expression levels of CCL3 in PBMCs
3. **Monocytes**: Cell frequency of monocytes in PBMCs

## Methodology

### Step 1: Data Harmonization

**Objective**: Standardize datasets from different years (2020-2023) for consistent analysis.

**Key Tasks**:
- Standardize data formats (units, feature names, sample identifiers)
- Feature matching across datasets
- Extract maximum overlapping information from datasets with varying feature coverage

### Step 2: Data Preprocessing

**Data Imputation**:
- Method: `impute.knn()` for handling missing data
- Rationale: Minimizes bias introduction during dataset merging

**Normalization**:
- Method: Median baseline value normalization
- Purpose: Account for measurement scale differences across datasets

**Batch Effect Correction**:
- Detection: Principal Component Analysis (PCA) for visualization
- Correction: ComBat algorithm implementation

### Step 3: Sample Availability Analysis

**Cross-Assay Sample Mapping**:
- Ensure consistent sample representation across all datasets
- Critical for multi-modal feature integration

**Assay Selection Strategy**:
- **Included**: Plasma Titer, PBMC Gene Expression, PBMC Cell Frequency
- **Excluded**: Cytokines (Olink/Legendplex), T Cell Activation/Polarization
- **Rationale**: Limited sample availability would require extensive imputation, potentially introducing bias and losing biological relevance

### Step 4: Single-Modal Predictive Modeling

**Approach**: Build baseline models using individual assay data.

**Model-Target Pairings**:
- IgG PT Titers ← Plasma Titer data
- CCL3 Gene Expression ← PBMC Gene Expression data  
- Monocyte Cell Frequency ← PBMC Cell Frequency data

**Modeling Pipeline**:
1. **Data Splitting**: Train/test split for performance evaluation
2. **Model Selection**: PyCaret or LazyPredict for automated model screening
3. **Model Training**: Train best-performing model on training set
4. **Model Evaluation**: Performance assessment using appropriate metrics
5. **Model Interpretation**: Feature importance and relationship analysis
6. **Model Tuning**: Hyperparameter optimization if needed

**Evaluation Metrics**:
- Accuracy, Precision, Recall, F1-score (classification)
- MAE, MSE, R² (regression)

### Step 5: Multi-Modal Predictive Modeling

**Trigger Condition**: Proceed if single-modal models show insufficient performance.

**Data Integration Strategy**:
- **Imputation Methods**: 
  - KNN imputation (fast, effective)
  - MICE (Multiple Imputation by Chained Equations) for complex structures
- **Integration Method**: MCIA (Multiple Co-Inertia Analysis)
  - Preserves inter-modal relationships
  - Identifies common patterns across modalities
  - Enhances predictive performance

**Pipeline**: Same as single-modal approach with integrated dataset input.
