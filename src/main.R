source("config.R")
source("src/import_data.R")
source("src/t_test.R")
source("src/anova.R")
source("src/metabolite_outliers.R")
source("src/metabolite_Normality_and_Equal_Variance_test.R")
source("src/cliff_delta.R")
source("src/multiple_linear_regression_mediation_analysis")

# 01-1. [EV] Check Conditions for T-Test
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

# 01-2. [EV] Perform Statistical Tests
#####

# Run T-tests
t_test_results_A <- t_test_report(input.BC.df[,2:4], "A")
t_test_results_B <- t_test_report(input.BC.df[,2:4], "B")

t_test_report(input.D.ex[,2:3], "C")
t_test_report(input.D.mi[,2:3], "C")

# Run ANOVA
anova_results_C <- anova_report(input.D.df[,2:3], "C")

# Save Results
write.csv(t_test_results_A, paste0(results_path, "t_test_results_A.csv"), row.names = FALSE)
write.csv(t_test_results_B, paste0(results_path, "t_test_results_B.csv"), row.names = FALSE)
write.csv(anova_results_C, paste0(results_path, "anova_results_C.csv"), row.names = FALSE)

# # Perform multiple statistical tests using LMSstat package
# library(LMSstat)
# 
# input.BC.stat <- Allstats(input.BC.df, Adjust_p_value = FALSE)
# input.BC.stat$Result
# input.D.stat <- Allstats(input.D.df, Adjust_p_value = FALSE)
# 
# # Extract statistical results
# input.D.stat$t_test[1,]
# input.D.stat$Anova
# input.D.stat$Anova_PostHoc[1,]
# 
# # Subset data and perform further statistical tests
# input.D.ex.stat <- Allstats(input.D.ex, Adjust_p_value = FALSE)
# input.D.ex.stat$t_test
# input.D.mi.stat <- Allstats(input.D.mi, Adjust_p_value = FALSE)
# input.D.mi.stat$t_test

#####


# 02-1. [Metabolomics] outlier test
#####
# Detect Outliers

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

# 02-2. [Metabolomics] Mann whitney u-test
#####
# metabolite normality, equal variance test
library(dplyr)
library(tidyr)
library(xlsx)
library(lawstat)

cell.output <- analyze_metabolites(cell.df)
media.output <- analyze_metabolites(media.df)
EV.output <- analyze_metabolites(EV.df)


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
#####

# 02-3. [Metabolomics] Cliff_delta cacluation
#####
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

# 02-4. [Metabolomics] multiple linear regression
#####
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

# 02-5. [Metabolomics] mediation analysis
#####
library(unglue)

nm.tmp <- asso_results.Group_Cell_Media.sig %>% filter(Pval <= 0.05)
nm <-  unique(c(unlist(nm.tmp["X"]), unlist(nm.tmp["Y"])))

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

final.results.01 <- mediation.results.tmp[,c(23,4:6,24:29,7:22)]

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

final.results.02 <- mediation.results.tmp[,c(23,4:6,24:29,7:22)]
                       
#####
