library(readr)
library(plotly)
library(ggplot2)
library(DescTools)
library(InformationValue)
library(tidyverse)
library(reshape2)
library(biglm)
library(Hmisc)
library(corrplot)
library(ggpubr)

#Reading the prs and phenotype data onto dataframes
dataframe <- function(SCZ_prs, BD_prs, WR_prs, EA_prs, DYX_prs, MDD_prs, covariates){
  #Reading the tsv files into dataframes
  SCZ_df = read.table(SCZ_prs, sep ="\t", header = TRUE)
  BD_df = read.table(BD_prs, sep="\t", header = TRUE)
  WR_df = read.table(WR_prs, sep="\t", header = TRUE)
  EA_df = read.table(EA_prs, sep="\t", header = TRUE)
  DYX_df = read.table(DYX_prs, sep="\t", header = TRUE)
  MDD_df = read.table(MDD_prs, sep="\t", header = TRUE)
  cov_df = read.table(covariates, sep="\t", header = FALSE)
  colnames(cov_df)=c("Participant_ID", "Age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Genotype_batch", "UK_Biobank_assessment_centre", "Date")
  cov_df[ ,c("UK_Biobank_assessment_centre", "Date")] <- list(NULL)
  
  #Renaming columns for ease of understanding
  colnames(SCZ_df) <- c('Participant_ID', 'FID', 'Allele_CT', 'Allele_Dosage_Sum', 'SCZ_PGS_Score')
  colnames(BD_df) <- c('Participant_ID', 'FID', 'Allele_CT', 'Allele_Dosage_Sum', 'BD_PGS_Score')
  colnames(WR_df) <- c('Participant_ID', 'FID', 'Allele_CT', 'Allele_Dosage_Sum', 'WR_PGS_Score')
  colnames(EA_df) <- c('Participant_ID', 'FID', 'Allele_CT', 'Allele_Dosage_Sum', 'EA_PGS_Score')
  colnames(DYX_df) <- c('Participant_ID', 'FID', 'Allele_CT', 'Allele_Dosage_Sum', 'DYX_PGS_Score')
  colnames(MDD_df) <- c('Participant_ID', 'FID', 'Allele_CT', 'Allele_Dosage_Sum', 'MDD_PGS_Score')
  
  return(list(SCZ = SCZ_df, BD = BD_df, WR = WR_df, EA = EA_df, DYX = DYX_df, MDD = MDD_df, cov = cov_df))
}

# Subset dataframe using common participants in all the dfs
make_merged_df <- function(dfs){
  
  # Only keeping the necessary columns
  dfs$SCZ = dfs$SCZ[,c('Participant_ID', 'SCZ_PGS_Score')]
  dfs$BD = dfs$BD[,c('Participant_ID', 'BD_PGS_Score')]
  dfs$WR = dfs$WR[,c('Participant_ID', 'WR_PGS_Score')]
  dfs$EA = dfs$EA[,c('Participant_ID', 'EA_PGS_Score')]
  dfs$DYX = dfs$DYX[,c('Participant_ID', 'DYX_PGS_Score')]
  dfs$MDD = dfs$MDD[,c('Participant_ID', 'MDD_PGS_Score')]
  
  # Changing character column to integer
  dfs$cov$Participant_ID = as.integer(dfs$cov$Participant_ID)
  
  # Merging columns based on Participation_ID
  merged_prs_df_w_cov <- dfs %>% reduce(full_join, by='Participant_ID')
  
  return(merged_prs_df_w_cov)
}

scaling <- function(prs_pheno_df){
  #Using standard scaling(mean value of zero, and SD = 1)
  prs_pheno_df$SCZ_PGS_Score <- scale(prs_pheno_df$SCZ_PGS_Score)
  prs_pheno_df$BD_PGS_Score <- scale(prs_pheno_df$BD_PGS_Score)
  prs_pheno_df$WR_PGS_Score <- scale(prs_pheno_df$WR_PGS_Score)
  prs_pheno_df$EA_PGS_Score <- scale(prs_pheno_df$EA_PGS_Score)
  prs_pheno_df$DYX_PGS_Score <- scale(prs_pheno_df$DYX_PGS_Score)
  prs_pheno_df$MDD_PGS_Score <- scale(prs_pheno_df$MDD_PGS_Score)
  
  # #Plotting the distribution of scores
  # density_fig <- ggplot(prs_pheno_df, aes(x=as.numeric(unlist(PRS_Score)), fill = as.factor(Qualification))) + 
  #   geom_density(alpha=0.4) +
  #   scale_fill_hue(labels = c("Low", "High")) +
  #   ggtitle("Distribution") +
  #   labs(fill = "Qualitfication", , x = "PRS Score", y = "Density")
  # density_fig
  # 
  # #Plotting a violin plot
  # violin_fig <- ggplot(prs_pheno_df, aes(x=as.factor(Qualification), y = as.numeric(unlist(PRS_Score)), , fill = as.factor(Qualification))) +
  #   geom_violin(alpha = 0.7) +
  #   geom_boxplot(width=0.2, color = "black") +
  #   labs(fill = "Qualitfication Status", x = "Qualification", y = "PRS Score") +
  #   scale_fill_hue(labels = c("Low", "High")) +
  #   ggtitle("Violon Plot")
  # violin_fig
  
  return(prs_pheno_df)
}

visual_graphs <- function(scaled_prs_df){
  
  subset <- scaled_prs_df[, c("SCZ_PGS_Score", "BD_PGS_Score", "WR_PGS_Score", "EA_PGS_Score")]
  data <- melt(subset)
  
  ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)
  
  plot1 <- ggplot(scaled_prs_df, aes(x=SCZ_PGS_Score))+
    geom_density(color="darkblue", fill="lightblue")
  plot1
  plot2 <- ggplot(scaled_prs_df, aes(x=BD_PGS_Score))+
    geom_density(color="magenta", fill="pink")
  plot2
  plot3 <- ggplot(scaled_prs_df, aes(x=WR_PGS_Score))+
    geom_density(color="green", fill="lightgreen")
  plot3
  plot4 <- ggplot(scaled_prs_df, aes(x=EA_PGS_Score))+
    geom_density(color="orange", fill="yellow")
  plot4
  
  plot5 <- ggplot(scaled_prs_df, aes(sample=SCZ_PGS_Score)) +
    stat_qq(color='blue') + 
    stat_qq_line()
  plot5
  plot6 <- ggplot(scaled_prs_df, aes(sample=BD_PGS_Score)) +
    stat_qq(color='red') +
    stat_qq_line()
  plot6
  plot7 <- ggplot(scaled_prs_df, aes(sample=WR_PGS_Score)) +
    stat_qq(color='green') +
    stat_qq_line()
  plot7
  plot8 <- ggplot(scaled_prs_df, aes(sample=EA_PGS_Score)) +
    stat_qq(color='orange') +
    stat_qq_line()
  plot8

  # plot.new()
  # tiff(file="distribution_new_1.tiff")
  ggarrange(plot1, plot2, plot3, plot4,
            labels = c("SCZ", "BD", "GenLang", "EA"),
            ncol = 4, nrow = 1)
  # dev.off()

  # ggarrange(plot5, plot6, plot7, plot8,
  #           labels = c("SCZ", "BD", "WR", "EA"),
  #           ncol = 4, nrow = 1)
}

correlation <- function(scaled_prs_df){
  
  scaled_prs_df <- na.omit(scaled_prs_df)
  
  cor_df <- scaled_prs_df[, c(2,3,4,5,6,7)]

  # result <- rcorr(as.matrix(cor_df), type=c("pearson", "spearman"))
  result <- Hmisc::rcorr(as.matrix(cor_df), type = c("pearson", "spearman"))

  # Extract correlation matrix and p-value matrix
  correlation_matrix <- result$r
  
  # Modify the row and column names (labels)
  new_labels <- c("SCZ", "BD", "WR", "EA", "DYX", "MDD")
  rownames(correlation_matrix) <- new_labels
  colnames(correlation_matrix) <- new_labels
  
  diag(result$P) <- 0
  p_value_matrix <- result$P
  
  # # Generate correlation plot
  # tiff(file = "corr_new_fixed.tiff", width = 12, height = 12, units = 'in', res = 150)
  # par(mar = c(4, 4, 4, 4))  # Adjust margins
  # 
  # # Main correlation plot with circles
  # corrplot::corrplot(
  #   correlation_matrix,
  #   method = "circle",        # Use circles
  #   type = "upper",           # Show upper triangle
  #   order = "hclust",         # Cluster similar correlations
  #   tl.col = "black",         # Label color
  #   tl.cex = 1.2,             # Text label size
  #   tl.srt = 45,              # Rotate text
  #   p.mat = p_value_matrix,   # Include p-values
  #   sig.level = 0.05,         # Only show significant correlations
  #   insig = "blank",          # Blank insignificant cells
  #   addCoef.col = "black",    # Add numeric correlation values
  #   col = colorRampPalette(c("red", "white", "blue"))(200)
  # )
  # dev.off()

  # plot.new()
  # tiff(file="corr_new.tiff", width = 12, height = 6, units = 'in', res = 150)
  # corrplot(result$r, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)
  # dev.off()
  return(result)
}

regression <- function(scaled_prs_df){
  model = lm(WR_PRS_Score~BD_PRS_Score + Age + sex, data = scaled_prs_df) #Create the linear regression
  #fit_model2 = lm(WR_PRS_Score~SCZ_PRS_Score + Age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Genotype_batch, data = scaled_prs_df) #Create the linear regression
  #return(list(m1 = fit_model1, m2 = fit_model2))
  #print(NagelkerkeR2(model))
  return(model)
}

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  plot.new()
  this_plot <- ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
  this_plot
}

SCZ_prs <- 'C:\\Project\\BP_SCZ_score\\SCZ_phi_auto.sscore'
BD_prs <- 'C:\\Project\\BP_SCZ_score\\BP_phi_auto.sscore'
WR_prs <- 'C:\\Project\\GenLang\\Genlang_auto.sscore'
EA_prs <- 'C:\\Project\\GenLang\\EA_auto.sscore'
DYX_prs <- 'C:\\Project\\DYX\\DYX_auto.sscore'
MDD_prs <- 'C:\\Project\\MDD\\MDD_auto.sscore'
covariates <- 'C:\\Project\\all_participants.tsv'

dfs <- dataframe(SCZ_prs, BD_prs, WR_prs, EA_prs, DYX_prs, MDD_prs, covariates)
prs_pheno_df <- make_merged_df(dfs)
scaled_prs_df <- scaling(prs_pheno_df)
visual_graphs(scaled_prs_df)
cor_res <- correlation(scaled_prs_df)
predicted_obs <- regression(scaled_prs_df)
ggplotRegression(predicted_obs)

tmp_df <- scaled_prs_df[, c(1,2,3,4,5)]
