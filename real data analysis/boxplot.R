#run data_preprocessing.R
# 
detect_outliers_iqr <- function(data, coef = 1.5) {
  outliers_list <- list()
  
  for(gene in colnames(data)) {
    values <- target_data_X[,which(colnames(target_data_X)==gene)]
    q <- quantile(values, probs = c(0.25, 0.75))
    iqr <- q[2] - q[1]
    lower_bound <- q[1] - coef * iqr
    upper_bound <- q[2] + coef * iqr
    outlier_indices <- which(values < lower_bound | values > upper_bound)
    outliers_list[[gene]] <- outlier_indices
  }
  
  return(outliers_list)
}
target_data_X = target_data$X
# 
outliers_iqr <- detect_outliers_iqr(target_data_X)
outlier_counts <- sapply(outliers_iqr, length)
top10_genes <- names(sort(outlier_counts, decreasing = TRUE))[1:10]

# 
top10_outlier_genes <- as.data.frame(target_data_X[, top10_genes])


library(tidyr)

# 
gd_long <- pivot_longer(top10_outlier_genes, cols = everything(), names_to = "gene", values_to = "expression")
# 
library(ggplot2)


# 
ggplot(data = gd_long, aes(x = gene, y = expression, fill = gene)) +
  geom_boxplot() + 
  labs(title = "Boxplot of Gene Expression", 
       x = "Gene", y = "Expression Level") +
  theme_minimal() 