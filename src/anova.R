anova_report <- function(data, variable) {
  aov_model <- aov(data[[variable]] ~ data$Group)
  anova_summary <- summary(aov_model)
  
  result <- data.frame(
    Variable = variable,
    DF_Between = anova_summary[[1]]$Df[1],
    DF_Within = anova_summary[[1]]$Df[2],
    F_Value = anova_summary[[1]]$`F value`[1],
    P_Value = anova_summary[[1]]$`Pr(>F)`[1]
  )
  return(result)
}