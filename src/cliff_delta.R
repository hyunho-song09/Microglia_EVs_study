library(effsize)

cliff_delta_fun <- function(data, variable) {
  group_0 <- data[[variable]][data$Group == "G0"]
  group_1 <- data[[variable]][data$Group == "G1"]
  
  delta_result <- cliff.delta(group_1, group_0)
  
  return(data.frame(
    Metric = variable,
    Delta = delta_result$estimate,
    Magnitude = as.character(delta_result$magnitude),
    P_Value = delta_result$p.value
  ))
}