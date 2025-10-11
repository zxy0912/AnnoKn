
########## check the risk region for different CHRs

s = numeric()

for(chrid in 1:22){
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/SCZ/annot_pvalus/risk_region/","EUR","_chr_", chrid,".txt")
  risk_region <- read.table(path)
  if(any(is.na(risk_region))){
    s = append(s, 0)
    next
  }
  print(dim(risk_region))
  s = append(s, nrow(risk_region))
}

chr <- 1:length(s)
df <- data.frame(Chromosome = factor(chr), Count = s)

library(ggplot2)

ggplot(df, aes(x = Chromosome, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Number of Risk Regions per Chromosome",
    x = "Chromosome",
    y = "Risk Region Count"
  )

########################

seed = '10'

region_1 <- numeric()
region_anno_1 <- numeric()
region_2 <- numeric()
region_anno_2 <- numeric()
lambdas1 <- list()
lambdas2 <- list()

threshold = 0.1
sum1 = 0
sum2 = 0
chrid = '1'
M = 1
# seed = '12345' for dss


chrlist = (1:22)[c(-21)]
for(chrid in chrlist){
  print(paste0("chrid", chrid))
  ancestry0 = 'EUR'
  ancestry1 = c('EAS',"AFR","LAT")
  # ancestry1 = c('EAS')
  other = paste(ancestry1, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/SCZ/annot_pvalus/result/", ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)
  
  print(length(result$risk_region))
  
  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)
  
  
  
  for(i in 1:n_region){
    # print(i)
    a = q_gk[[i]]
    if(sum(a <= threshold) > 0){
      region_1 <- rbind(region_1, as.numeric(append(chrid, result$risk_region[[i]])))
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      region_anno_1 <- rbind(region_anno_1, as.numeric(append(chrid, result$risk_region[[i]])))
      lambdas1 <- append(lambdas1, list(result$lambda_s[[i]]))
    }
  }
  
  
  ancestry0 = 'EAS'
  ancestry1 = c('EUR',"AFR","LAT")
  # ancestry1 = c('EUR')
  other = paste(ancestry1, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/SCZ/annot_pvalus/result/", ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)
  print(length(result$risk_region))
  
  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)
  
  
  for(i in 1:n_region){
    # print(i)
    a = q_gk[[i]]
    if(sum(a <= threshold) > 0){
      region_2 <- rbind(region_2, as.numeric(append(chrid, result$risk_region[[i]])))
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      region_anno_2 <- rbind(region_anno_2, as.numeric(append(chrid, result$risk_region[[i]])))
      lambdas2 <- append(lambdas2, list(result$lambda_s[[i]]))
    }
  }
  
}

nrow(region_1)
nrow(region_2)

nrow(region_anno_1)
nrow(region_anno_2)




df1 <- as.data.frame(region_1)
df2 <- as.data.frame(region_2)
common_rows <- merge(df1, df2)
n_common <- nrow(common_rows)
n_common


df1 <- as.data.frame(region_anno_1)
df2 <- as.data.frame(region_anno_2)
common_rows <- merge(df1, df2)
n_common <- nrow(common_rows)
n_common


lambda_matrix1 <- do.call(rbind, lambdas1)
lambda_matrix2 <- do.call(rbind, lambdas2)

j = 2

library(ggplot2)
df <- data.frame(lambdas = lambda_matrix1[,j])
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,  
    fill = "#0072B2", 
    color = "white",
    alpha = 0.8     
  ) +
  labs(
    title = "Distribution of lambda_s of AFR p-values for EUR GWAS",
    x = "Lambda_s",            
    y = "Frequency (Count)"         
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot



df <- data.frame(lambdas = lambda_matrix2[,j])
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,  
    fill = "#0072B2", 
    color = "white",
    alpha = 0.8     
  ) +
  labs(
    title = "Distribution of lambda_s of AFR p-values for EAS GWAS",
    x = "Lambda_s",            
    y = "Frequency (Count)"         
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot



############## check if multiple knockoff copies


find1 <- c()
find1_anno <- c()
find2 <- c()
find2_anno <- c()


for(seed in c(1,12,10,10000,10000, 100, 1000, 123, 1234, 12345)){
  region_1 <- numeric()
  region_anno_1 <- numeric()
  region_2 <- numeric()
  region_anno_2 <- numeric()
  lambdas1 <- list()
  lambdas2 <- list()
  
  threshold = 0.1
  sum1 = 0
  sum2 = 0
  chrid = '1'
  M = 1
  # seed = '12345' for dss
  
  
  chrlist = (1:22)[c(-21)]
  for(chrid in chrlist){
    print(paste0("chrid", chrid))
    ancestry0 = 'EUR'
    ancestry1 = c('EAS',"AFR","LAT")
    # ancestry1 = c('EAS')
    other = paste(ancestry1, collapse = "_")
    path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/SCZ/annot_pvalus/result/", ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    load(path)
    
    print(length(result$risk_region))
    
    q_gk = result$q_gk
    q_gkanno = result$q_gkanno
    n_region = length(q_gk)
    
    
    
    for(i in 1:n_region){
      # print(i)
      a = q_gk[[i]]
      if(sum(a <= threshold) > 0){
        region_1 <- rbind(region_1, as.numeric(append(chrid, result$risk_region[[i]])))
        find1 <- append(find1, paste(chrid, i, sep = '_'))
      }
      b = q_gkanno[[i]]
      # print(which(b <= threshold))
      if(sum(b <= threshold) > 0){
        region_anno_1 <- rbind(region_anno_1, as.numeric(append(chrid, result$risk_region[[i]])))
        lambdas1 <- append(lambdas1, list(result$lambda_s[[i]]))
        find1_anno <- append(find1_anno, paste(chrid, i, sep = '_'))
      }
    }
    
    
    ancestry0 = 'EAS'
    ancestry1 = c('EUR',"AFR","LAT")
    # ancestry1 = c('EUR')
    other = paste(ancestry1, collapse = "_")
    path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/SCZ/annot_pvalus/result/", ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
    load(path)
    print(length(result$risk_region))
    
    q_gk = result$q_gk
    q_gkanno = result$q_gkanno
    n_region = length(q_gk)
    
    
    for(i in 1:n_region){
      # print(i)
      a = q_gk[[i]]
      if(sum(a <= threshold) > 0){
        region_2 <- rbind(region_2, as.numeric(append(chrid, result$risk_region[[i]])))
        find2 <- append(find2, paste(chrid, i, sep = '_'))
      }
      b = q_gkanno[[i]]
      # print(which(b <= threshold))
      if(sum(b <= threshold) > 0){
        region_anno_2 <- rbind(region_anno_2, as.numeric(append(chrid, result$risk_region[[i]])))
        lambdas2 <- append(lambdas2, list(result$lambda_s[[i]]))
        find2_anno <- append(find2_anno, paste(chrid, i, sep = '_'))
      }
    }
    
  }
  
  nrow(region_1)
  nrow(region_2)
  
  nrow(region_anno_1)

}

repli = 4

sum(table(find1) > repli)
sum(table(find1_anno) > repli)

sum(table(find2) > repli)
sum(table(find2_anno) > repli)


length(table(find1))
length(table(find1_anno))
length(table(find2))
length(table(find2_anno))



length(intersect(unique(find1), unique(find2)))
length(intersect(unique(find1_anno), unique(find2_anno)))

unique(find1)[order(unique(find1))]
unique(find2)[order(unique(find2))]

unique(find1_anno)[order(unique(find1_anno))]
unique(find2_anno)[order(unique(find2_anno))]















library(ggplot2)
library(dplyr)
library(tidyr)

ghost_counts <- as.numeric(table(find2))
anno_counts  <- as.numeric(table(find2_anno))

thresholds <- 1:10

ghost_summary <- data.frame(
  threshold = thresholds/10,
  count = sapply(thresholds, function(t) sum(ghost_counts >= t)),
  Method = "GhostKnockoff"
)

anno_summary <- data.frame(
  threshold = thresholds/10,
  count = sapply(thresholds, function(t) sum(anno_counts >= t)),
  Method = "AnnoGK"
)

plot_df <- rbind(ghost_summary, anno_summary)

p_thresh <- ggplot(plot_df, aes(x = factor(threshold), y = count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Number of signals detected ≥ threshold (EAS GWAS)",
    x = "Detection frequency threshold",
    y = "Number of signals"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank()
  )

p_thresh
