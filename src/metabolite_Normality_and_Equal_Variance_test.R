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
    
    # levene.test
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