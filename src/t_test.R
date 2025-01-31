library(dplyr)

t_test_report <- function(data, variable) {
  t_test <- t.test(data[[variable]] ~ data$Group, var.equal = TRUE)
  
  result <- data.frame(
    Variable = variable,
    DF = t_test$parameter,
    T_Value = t_test$statistic,
    P_Value = t_test$p.value
  )
  return(result)
}