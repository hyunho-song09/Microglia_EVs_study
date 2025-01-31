detect_outliers_robust <- function(data, group_col, group, robust_thresh = 5) {
  group_data <- data[data[[group_col]] == group, ]
  
  outlier_counts <- data.frame(metabolite = colnames(group_data)[3:ncol(group_data)], outlier_count = 0)
  
  for (met in colnames(group_data)[3:ncol(group_data)]) {
    values <- group_data[[met]]
    median_val <- median(values, na.rm = TRUE)
    sd_val <- sd(values, na.rm = TRUE)
    
    robust_scores <- (values - median_val) / sd_val
    
    outlier_counts$outlier_count[outlier_counts$metabolite == met] <- sum(abs(robust_scores) > robust_thresh)
  }
  return(outlier_counts)
}