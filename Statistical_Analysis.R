################################
########## EV Part #############
################################


# 01-1. [EV] Import Data
#####
rm(list=ls())

library(xlsx)

# Load datasets
input.BC.df <- read.xlsx("example_data_01.xlsx", sheetName = "input_fig_BC")
input.D.df <- read.xlsx("example_data_01.xlsx", sheetName = "input_fig_D")
input.D.ex <- input.D.df[c(1:3,7:9),]
input.D.mi <- input.D.df[c(4:6,10:12),]
#####

# 01-2. [EV] Check Conditions for T-Test
#####
library(lawstat)

# Normality test for each group
by(input.BC.df$A, input.BC.df$Group, shapiro.test)
by(input.BC.df$B, input.BC.df$Group, shapiro.test)
by(input.D.df$C, input.D.df$Group, shapiro.test)

# Homogeneity of variance tests

# Figure B
levene.test(input.BC.df$A, input.BC.df$Group)$p.value
bartlett.test(input.BC.df$A, input.BC.df$Group)$p.value

# Figure C
levene.test(input.BC.df$B, input.BC.df$Group)$p.value
bartlett.test(input.BC.df$B, input.BC.df$Group)$p.value

# Figure D
levene.test(input.D.df$C, input.D.df$Group)$p.value
bartlett.test(input.D.df$C, input.D.df$Group)$p.value

# Additional variance tests for Figure D ~ F
levene.test(input.D.df[c(1:3,7:9),]$C, input.D.df[c(1:3,7:9),]$Group)$p.value
levene.test(input.D.df[c(4:6,10:12),]$C, input.D.df[c(4:6,10:12),]$Group)$p.value
bartlett.test(input.D.df[c(1:3,7:9),]$C, input.D.df[c(1:3,7:9),]$Group)$p.value
bartlett.test(input.D.df[c(4:6,10:12),]$C, input.D.df[c(4:6,10:12),]$Group)$p.value
#####

# 01-3. [EV] Perform Statistical Tests
#####

# Load necessary library
library(dplyr)

# Function for performing t-test and reporting results
t_test_report <- function(data, variable) {
  t_test <- t.test(data[[variable]] ~ data$Group, var.equal = TRUE)
  
  cat("\nT-test for", variable, "\n")
  cat("Degrees of Freedom:", t_test$parameter, "\n")
  cat("T value:", t_test$statistic, "\n")
  cat("P value:", t_test$p.value, "\n")
}

# Perform t-tests
t_test_report(input.BC.df[,2:4], "A")
t_test_report(input.BC.df[,2:4], "B")

t_test_report(input.D.ex[,2:3], "C")
t_test_report(input.D.mi[,2:3], "C")

# Function for performing ANOVA and reporting results
anova_report <- function(data, variable) {
  aov_model <- aov(data[[variable]] ~ data$Group)
  anova_summary <- summary(aov_model)
  
  cat("\nANOVA for", variable, "\n")
  cat("Degrees of Freedom: Between =", anova_summary[[1]]$Df[1], ", Within =", anova_summary[[1]]$Df[2], "\n")
  cat("F value:", anova_summary[[1]]$`F value`[1], "\n")
  cat("P value:", anova_summary[[1]]$`Pr(>F)`[1], "\n")
}

# Perform ANOVA
anova_report(input.D.df[,2:3], "C")

# Perform multiple statistical tests using LMSstat package
library(LMSstat)

input.BC.stat <- Allstats(input.BC.df, Adjust_p_value = FALSE)
input.BC.stat$Result

input.D.stat <- Allstats(input.D.df, Adjust_p_value = FALSE)

# Extract statistical results
input.D.stat$t_test[1,]
input.D.stat$Anova
input.D.stat$Anova_PostHoc[1,]

# Subset data and perform further statistical tests
input.D.ex.stat <- Allstats(input.D.ex, Adjust_p_value = FALSE)
input.D.ex.stat$t_test
input.D.mi.stat <- Allstats(input.D.mi, Adjust_p_value = FALSE)
input.D.mi.stat$t_test

#####


################################
###### Metabolomics Part #######
################################

# 02-1. [Metabolomics] Import Data
#####
rm(list=ls())

setwd("D:/project/experiment/BrainOmics/AD_omics/241024_EV_manuscript")

library(xlsx)
cell.df <- read.xlsx("example_data_02.xlsx", sheetName = "ex2")
media.df <- read.xlsx("example_data_03.xlsx", sheetName = "ex3")
EV.df <- read.xlsx("example_data_04.xlsx", sheetName = "ex4")
#####

# 02-2. [Metabolomics] outlier test
#####
# Perform Robust Standardization based outlier detection
detect_outliers_robust <- function(data, group_col, group, robust_thresh = 5) {
  # Subset data by group
  group_data <- data[data[[group_col]] == group, ]
  
  # Initialize outlier count dataframe
  outlier_counts <- data.frame(metabolite = colnames(group_data)[3:ncol(group_data)], outlier_count = 0)
  
  # Calculate Robust Standardization scores for each metabolite
  for (met in colnames(group_data)[3:ncol(group_data)]) {
    values <- group_data[[met]]
    median_val <- median(values, na.rm = TRUE) # Compute median
    sd_val <- sd(values, na.rm = TRUE)         # Compute standard deviation
    
    # Robust standardization
    robust_scores <- (values - median_val) / sd_val
    
    # Count outliers based on threshold
    outlier_counts$outlier_count[outlier_counts$metabolite == met] <- sum(abs(robust_scores) > robust_thresh)
  }
  
  return(outlier_counts)
}

# Apply Robust Standardization to different datasets

# cell
ab_outliers_cell <- detect_outliers_robust(cell.df, group_col = "Group", group = "G1", robust_thresh = 3)
control_outliers_cell <- detect_outliers_robust(cell.df, group_col = "Group", group = "G0", robust_thresh = 3)
cell_outlier <- merge(ab_outliers_cell, control_outliers_cell, by = "metabolite", suffixes = c("_G1", "_G0"))

# media
ab_outliers_media <- detect_outliers_robust(media.df, group_col = "Group", group = "G1", robust_thresh = 3)
control_outliers_media <- detect_outliers_robust(media.df, group_col = "Group", group = "G0", robust_thresh = 3)
media_outlier <- merge(ab_outliers_media, control_outliers_media, by = "metabolite", suffixes = c("_G1", "_G0"))

# EV
ab_outliers_EV <- detect_outliers_robust(EV.df, group_col = "Group", group = "G1", robust_thresh = 3)
control_outliers_EV <- detect_outliers_robust(EV.df, group_col = "Group", group = "G0", robust_thresh = 3)
EV_outlier <- merge(ab_outliers_EV, control_outliers_EV, by = "metabolite", suffixes = c("_G1", "_G0"))
#####

# 02-3. [Metabolomics] Mann whitney u-test
#####

library(LMSstat)

cell.stats <- Allstats(cell.df, Adjust_p_value = F)
cell.results <- data.frame(cell.stats$Result,
                           adj_ttest = p.adjust(cell.stats$t_test, method = "BH"),
                           adj_utest = p.adjust(cell.stats$u_test, method = "BH"))

media.stats <- Allstats(media.df, Adjust_p_value = F)
media.results <- data.frame(media.stats$Result,
                            p.adjust(media.stats$t_test, method = "BH"),
                            p.adjust(media.stats$u_test, method = "BH"))


EV.stats <- Allstats(EV.df, Adjust_p_value = F)
EV.results <- data.frame(EV.stats$Result,
                         p.adjust(EV.stats$t_test, method = "BH"),
                         p.adjust(EV.stats$u_test, method = "BH"))

# metabolite normality, equal variance test
library(dplyr)
library(tidyr)
library(xlsx)
library(lawstat)

analyze_metabolites <- function(data) {
  results <- data.frame(Metabolite = character(), 
                        Normality_AB_p = numeric(), 
                        Normality_Control_p = numeric(), 
                        Normality = integer(), 
                        Homogeneity_bartlett_p = numeric(), 
                        Homogeneity_bartlett = integer(),
                        Homogeneity_levene_p = numeric(), 
                        Homogeneity_levene = integer(), 
                        stringsAsFactors = FALSE)
  
  group_ab <- data %>% filter(Group == 'G1')
  group_control <- data %>% filter(Group == 'G0')
  
  for (met in colnames(data)[3:ncol(data)]) {
    ab_values <- group_ab[[met]]
    control_values <- group_control[[met]]
    
    shapiro_ab <- shapiro.test(ab_values)$p.value
    shapiro_control <- shapiro.test(control_values)$p.value
    
    normality <- ifelse(shapiro_ab > 0.05 & shapiro_control > 0.05, 1, 0)
    
    if (normality == 1) {
      bartlett_p <- bartlett.test(list(ab_values, control_values))$p.value
      homogeneity_bartlett <- ifelse(bartlett_p > 0.05, 1, 0)
    } else {
      bartlett_p <- NA
      homogeneity_bartlett <- 0
    }
    
    # levene.test 수정
    levene_p <- levene.test(c(ab_values, control_values), 
                            group = rep(1:2, times = c(length(ab_values), length(control_values))))$p.value
    
    homogeneity_levene <- ifelse(levene_p > 0.05, 1, 0)
    
    results <- rbind(results, data.frame(
      Metabolite = met, 
      Normality_AB_p = shapiro_ab, 
      Normality_Control_p = shapiro_control, 
      Normality = normality, 
      Homogeneity_bartlett_p = bartlett_p, 
      Homogeneity_bartlett = homogeneity_bartlett,
      Homogeneity_levene_p = levene_p, 
      Homogeneity_levene = homogeneity_levene
    ))
  }
  return(results)
}

cell.output <- analyze_metabolites(cell.df)
media.output <- analyze_metabolites(media.df)
EV.output <- analyze_metabolites(EV.df)

#####

# 02-4. [Metabolomics] Cliff_delta cacluation
#####

library(effsize)

cliff_delta_fun <- function(data, variable) {
  # Extract values for each group
  group_0 <- data[[variable]][data$Group == "G0"]
  group_1 <- data[[variable]][data$Group == "G1"]
  
  # Compute Cliff's Delta
  delta_result <- cliff.delta(group_1, group_0)
  
  # Return a data frame with results
  return(data.frame(
    Metric = variable,
    Delta = delta_result$estimate,
    Magnitude = as.character(delta_result$magnitude),
    stringsAsFactors = FALSE
  ))
}

cell.results <- data.frame()
for (col in colnames(cell.df)[3:ncol(cell.df)]){
  cell.results <- rbind(cell.results, cliff_delta_fun(cell.df, col))}

media.results <- data.frame()
for (col in colnames(media.df)[3:ncol(media.df)]){
  media.results <- rbind(media.results, cliff_delta_fun(media.df, col))}

EV.results <- data.frame()
for (col in colnames(EV.df)[3:ncol(EV.df)]){
  EV.results <- rbind(EV.results, cliff_delta_fun(EV.df, col))}

#####

# 02-5. [Metabolomics] Multiple linear regression
#####

library(dplyr)
library(doParallel)
library(doSNOW)

# creat function
meta.MLR.calculator <- function (X,Y,C,covariate) {
  print("Running...")
  
  # Internal Use Only / Private Code
  
}

# 01) Group vs Cell
MLR_Group_Cell.sig <- meta.MLR.calculator(X=Group.df.sig,
                                          Y=Cell.df.sig,
                                          C=NULL,
                                          covariate=FALSE)
MLR_Group_Cell.sig$p.adj_BH <- p.adjust(MLR_Group_Cell.sig$Pval, method="BH")
MLR_Group_Cell.sig$x_class <- "Group"
MLR_Group_Cell.sig$y_class <- "Cell"

# 02) Cell vs Media
MLR_Cell_Media.sig <- meta.MLR.calculator(X=Cell.df.sig,
                                          Y=Media.df.sig,
                                          C=NULL,
                                          covariate=FALSE)
MLR_Cell_Media.sig$p.adj_BH <- p.adjust(MLR_Cell_Media.sig$Pval, method="BH")
MLR_Cell_Media.sig$x_class <- "Cell"
MLR_Cell_Media.sig$y_class <- "Media"

# 03) Group vs Media
MLR_Group_Media.sig <- meta.MLR.calculator(X=Group.df.sig,
                                           Y=Media.df.sig,
                                           C=NULL,
                                           covariate=FALSE)
MLR_Group_Media.sig$p.adj_BH <- p.adjust(MLR_Group_Media.sig$Pval, method="BH")
MLR_Group_Media.sig$x_class <- "Group"
MLR_Group_Media.sig$y_class <- "Media"


# 04) Cell vs EV
MLR_Cell_EV.sig <- meta.MLR.calculator(X=Cell.df.sig,
                                       Y=EV.df.sig,
                                       C=NULL,
                                       covariate=FALSE)
MLR_Cell_EV.sig$p.adj_BH <- p.adjust(MLR_Cell_EV.sig$Pval, method="BH")
MLR_Cell_EV.sig$x_class <- "Cell"
MLR_Cell_EV.sig$y_class <- "EV"

# 05) Group vs EV
MLR_Group_EV.sig <- meta.MLR.calculator(X=Group.df.sig,
                                        Y=EV.df.sig,
                                        C=NULL,
                                        covariate=FALSE)
MLR_Group_EV.sig$p.adj_BH <- p.adjust(MLR_Group_EV.sig$Pval, method="BH")
MLR_Group_EV.sig$x_class <- "Group"
MLR_Group_EV.sig$y_class <- "EV"

# Group vs Cell vs Media
asso_results.tmp.sig <- rbind(MLR_Group_Cell.sig, MLR_Cell_Media.sig, MLR_Group_Media.sig)
asso_results.Group_Cell_Media.sig <- asso_results.tmp.sig %>% filter(Pval <= 0.05)
# Group vs Cell vs EV
asso_results.tmp.sig <- rbind(MLR_Group_Cell.sig, MLR_Cell_EV.sig, MLR_Group_EV.sig)
asso_results.Group_Cell_EV.sig <- asso_results.tmp.sig %>% filter(Pval <= 0.05)

#####

# 02-6. [Metabolomics] Mediation analysis
#####
library(unglue)

nm.tmp <- asso_results.Group_Cell_Media.sig %>% filter(Pval <= 0.05)
nm <-  unique(c(unlist(nm.tmp["X"]), unlist(nm.tmp["Y"])))

library(dplyr)
library(doParallel)
library(doSNOW)

# create function
meta.mediate.calculator <- function (X,M,Y,C,covariate) {
  print("Running...")
  
  # Internal Use Only / Private Code
  
}

# output media
# Create input dataframe
med.input.df.01 <- Group.df.sig
med.input.df.02 <- Cell.df.sig[,colnames(Cell.df.sig) %in% nm]
med.input.df.03 <- Media.df.sig[,colnames(Media.df.sig) %in% nm]

mediation.results.01 <- meta.mediate.calculator(X=med.input.df.01,
                                                M=med.input.df.02,
                                                Y=med.input.df.03,
                                                C=NULL,
                                                covariate=FALSE)

# integrate Multiple linear regression + mediation
asso_results.tmp <- asso_results.Group_Cell_Media.sig
asso_results.tmp$X_M_index <- paste0(asso_results.tmp$X, asso_results.tmp$Y)
asso_results.tmp$X_Y_index <- asso_results.tmp$X_M_index
asso_results.tmp$M_Y_index <- asso_results.tmp$X_M_index

mediation.results.tmp <- mediation.results.01
mediation.results.tmp$index <-  1:nrow(mediation.results.tmp)
mediation.results.tmp$X_M_index <- paste0(mediation.results.tmp$X, mediation.results.tmp$M)
mediation.results.tmp$X_Y_index <- paste0(mediation.results.tmp$X, mediation.results.tmp$Y)
mediation.results.tmp$M_Y_index <- paste0(mediation.results.tmp$M, mediation.results.tmp$Y)

mediation.results.tmp <- merge(mediation.results.tmp, asso_results.tmp[,c("X_M_index","coef","Pval")], by = "X_M_index")
colnames(mediation.results.tmp)[24:25] <- c("X_M_coef", "X_M_pval")
mediation.results.tmp <- merge(mediation.results.tmp, asso_results.tmp[,c("X_Y_index","coef","Pval")], by = "X_Y_index")
colnames(mediation.results.tmp)[26:27] <- c("X_Y_coef", "X_Y_pval")
mediation.results.tmp <- merge(mediation.results.tmp, asso_results.tmp[,c("M_Y_index","coef","Pval")], by = "M_Y_index")
colnames(mediation.results.tmp)[28:29] <- c("M_Y_coef", "M_Y_pval")

final.results <- cbind(mediation.results.tmp[,c(23,4:6,24:29,7:22)])


# output EV
# Create input dataframe
nm.tmp <- asso_results.Group_Cell_EV.sig %>% filter(Pval <= 0.05)
nm <-  unique(c(unlist(nm.tmp["X"]), unlist(nm.tmp["Y"])))

med.input.df.01 <- Group.df.sig
med.input.df.02 <- Cell.df.sig[,colnames(Cell.df.sig) %in% nm]
med.input.df.03 <- EV.df.sig[,colnames(EV.df.sig) %in% nm]

mediation.results.01 <- meta.mediate.calculator(X=med.input.df.01,
                                                M=med.input.df.02,
                                                Y=med.input.df.03,
                                                C=NULL,
                                                covariate=FALSE)

# integrate Multiple linear regression + mediation
asso_results.tmp <- asso_results.Group_Cell_EV.sig
asso_results.tmp$X_M_index <- paste0(asso_results.tmp$X, asso_results.tmp$Y)
asso_results.tmp$X_Y_index <- asso_results.tmp$X_M_index
asso_results.tmp$M_Y_index <- asso_results.tmp$X_M_index

mediation.results.tmp <- mediation.results.01
mediation.results.tmp$index <-  1:nrow(mediation.results.tmp)
mediation.results.tmp$X_M_index <- paste0(mediation.results.tmp$X, mediation.results.tmp$M)
mediation.results.tmp$X_Y_index <- paste0(mediation.results.tmp$X, mediation.results.tmp$Y)
mediation.results.tmp$M_Y_index <- paste0(mediation.results.tmp$M, mediation.results.tmp$Y)

mediation.results.tmp <- merge(mediation.results.tmp, asso_results.tmp[,c("X_M_index","coef","Pval")], by = "X_M_index")
colnames(mediation.results.tmp)[24:25] <- c("X_M_coef", "X_M_pval")
mediation.results.tmp <- merge(mediation.results.tmp, asso_results.tmp[,c("X_Y_index","coef","Pval")], by = "X_Y_index")
colnames(mediation.results.tmp)[26:27] <- c("X_Y_coef", "X_Y_pval")
mediation.results.tmp <- merge(mediation.results.tmp, asso_results.tmp[,c("M_Y_index","coef","Pval")], by = "M_Y_index")
colnames(mediation.results.tmp)[28:29] <- c("M_Y_coef", "M_Y_pval")

final.results <- cbind(mediation.results.tmp[,c(23,4:6,24:29,7:22)])


#####






