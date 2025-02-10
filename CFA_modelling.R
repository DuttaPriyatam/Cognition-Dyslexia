# Comprehensive analysis of cognitive factors
library(lavaan)
library(dplyr)
library(tidyr)
library(ggplot2)

# Function to extract and document factor scores with full transparency
extract_documented_scores <- function(data) {
  # Specify the model
  model <- '
    # # Speed factor
    # Speed =~ PM_spd + PRM_spd + SDS_spd + TM_num_spd + TM_anum_spd + RT_spd
    # 
    # # Accuracy factor
    # Accuracy =~ RT_acc + PM_acc + PRM_acc + FI_acc + NM_acc + SDS_acc
    # 
    # # # Residual correlations for same tasks
    # PM_spd ~~ PM_acc
    # PRM_spd ~~ PRM_acc
    # SDS_spd ~~ SDS_acc
    # 
    # # Additional correlations for similar tasks
    # TM_num_spd ~~ TM_anum_spd
    # FI_acc ~~ NM_acc
    # 
    # # Factor correlation
    # Speed ~~ Accuracy
    
    General =~ PM_spd + PRM_spd + SDS_spd + TM_num_spd + TM_anum_spd +
    RT_acc + PM_acc + PRM_acc + FI_acc + NM_acc + SDS_acc + RT_spd

    # Specific factors
    Speed =~ PM_spd + PRM_spd + SDS_spd + TM_num_spd + TM_anum_spd + RT_spd
    Accuracy =~ RT_acc + PM_acc + PRM_acc + FI_acc + NM_acc + SDS_acc

    # Orthogonality constraints for bifactor structure
    General ~~ 0*Speed
    General ~~ 0*Accuracy
    Speed ~~ 0*Accuracy

    # Residual correlations
    RT_spd ~~ RT_acc
    PM_spd ~~ PM_acc
    SDS_spd ~~ SDS_acc
    TM_num_spd ~~ TM_anum_spd
  '
  
  # Fit the model
  fit <- cfa(model,
             data = data,
             estimator = "MLR",
             missing = "fiml",
             std.lv = TRUE)
  
  # Extract fit measures
  fit_indices <- fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr", 
                                    "aic", "bic", "chisq", "df", "pvalue"))
  
  # Get standardized solution and R² values
  std_solution <- standardizedSolution(fit)
  loadings <- std_solution %>%
    filter(op == "=~") %>%
    mutate(
      R2 = est.std^2,
      Factor = lhs,
      Indicator = rhs
    ) %>%
    select(Factor, Indicator, est.std, R2, se, z, pvalue) %>%
    arrange(Factor, desc(R2))
  
  # Calculate factor scores
  scores <- lavPredict(fit, type = "lv")
  scores_df <- as.data.frame(scores)
  
  # Calculate factor score reliability
  r2_by_factor <- loadings %>%
    group_by(Factor) %>%
    summarise(
      Mean_R2 = mean(R2),
      Min_R2 = min(R2),
      Max_R2 = max(R2),
      Total_R2 = sum(R2),
      N_indicators = n()
    )
  
  # Return all results
  return(list(
    scores = scores_df,
    fit_indices = fit_indices,
    loadings = loadings,
    factor_reliability = r2_by_factor,
    fit = fit
  ))
}

# Example usage:
data <- read.table('cog_res.tsv', sep ="\t", header = TRUE)
results <- extract_documented_scores(data)

# Print comprehensive results
cat("MODEL FIT INDICES:\n")
print(results$fit_indices)

cat("\nFACTOR LOADINGS AND R² VALUES:\n")
print(results$loadings)

cat("\nFACTOR RELIABILITY SUMMARY:\n")
print(results$factor_reliability)

# Add ID column to scores if it doesn't exist
scores_with_id <- data.frame(
  Speed = results$scores$Speed,
  Accuracy = results$scores$Accuracy
)

# Merge with original data
merged_data <- cbind(data, 
                     Speed = results$scores$Speed,
                     Accuracy = results$scores$Accuracy,
                     General = results$scores$General)

# Save as TSV
write.table(merged_data,
            "cognitive_data_with_factors_final_final.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

# Print first few rows to verify
head(merged_data)

# # Print dimensions to verify
# dim(merged_data)


# Check correlation between Speed and Accuracy factors
cor(results$scores$Speed, results$scores$Accuracy, use = "pairwise.complete.obs")

# Extract latent variable covariance matrix
lavInspect(results$fit, "cov.lv")

# Get factor correlation
standardizedSolution(results$fit) %>%
  filter(op == "~~", lhs == "Speed", rhs == "Accuracy")
