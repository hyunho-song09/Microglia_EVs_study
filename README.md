# 🧬 Microglia EVs study

This repository contains **R scripts** and **sample data** used for metabolomics and extracellular vesicles (EVs) analysis in primary microglia cell related Alzheimer's Disease (AD) research. The pipeline includes data preprocessing, statistical analysis, multiple linear regression, and mediation analysis.
<br />


## 📖 Overview

This project implements:
- **Data Import & Preprocessing**: Reads Excel-based omics datasets.
- **Statistical Testing**: T-tests, ANOVA, normality, and variance tests.
- **Metabolomics Analysis**: Outlier detection, Cliff’s Delta calculation.
- **Multiple Linear Regression (MLR)**: Identifies relationships between variables.
- **Mediation Analysis**: Explores potential mediators in metabolic pathways.
<br />

---
## 🚀 How to Run
To execute the full analysis pipeline, run:
```
source("main.R")
```
---
<br />


## 📊 Pipeline Details
### 1️⃣ Data Import & Preprocessing
#### 📂 Scripts Used:
``` src/import_data.R ```
``` src/metabolite_outliers.R ```
#### Description:
- Loads metabolomics and EV data from .xlsx files.
- Identifies and removes outliers using robust standardization.
<br />

### 2️⃣ Statistical Analysis
#### 📂 Scripts Used:
``` src/t_test.R ```
``` src/anova.R ```
``` src/metabolite_Normality_and_Equal_Variance_test.R ```
#### Description:
- T-tests and ANOVA to compare groups.
- Shapiro-Wilk test for normality.
- Levene’s & Bartlett’s tests for homogeneity of variance.
<br />

### 3️⃣ Cliff’s Delta Calculation
#### 📂 Script Used:
``` src/cliff_delta.R ```
#### Description:
- Calculates Cliff’s Delta, a non-parametric effect size measure
<br />

### 4️⃣ Multiple Linear Regression (MLR)
#### 📂 Script Used: 
``` src/multiple_linear_regression_mediation_analysis.R ```
#### Description:
- Performs multiple linear regression (MLR) across datasets:
1) Group vs Cell
2) Cell vs Media
3) Group vs Media
4) Cell vs EV
5) Group vs EV
- Adjusts p-values using Benjamini-Hochberg correction.
<br />

### 5️⃣ Mediation Analysis
#### 📂 Script Used: 
``` src/multiple_linear_regression_mediation_analysis.R ```
#### Description:
- Identifies metabolic mediators between variables.
- Integrates MLR results with mediation models.

**Seong-Hun Bong, Hayoung Choi, Hyun-Ho Song, Dong Kyu Kim, Inhee Mook-Jung, Do Yup Lee** (2025).  
*Metabolic Reprogramming in Primary Microglial Cell and Extracellular Vesicle Triggered by Aβ Exposure*.  
**Journal of Neurochemistry**. https://doi.org/10.1111/jnc.70030
