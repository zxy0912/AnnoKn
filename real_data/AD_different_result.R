
########## check the risk region for different CHRs
s = numeric()

for(chrid in 1:22){
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/AD_new/annot_pvalus/risk_region/","GWAS1","_chr_", chrid,".txt")
  risk_region <- read.table(path)
  if(any(is.na(risk_region))){
    s = append(s, 0)
    next
  }
  print(dim(risk_region))
  s = append(s, nrow(risk_region))
}
########################

seed = '123'

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


chrlist = (1:22)
for(chrid in chrlist){
  
  if(s[chrid] == 0){
    next
  }
  print(paste0("chrid", chrid))
  ancestry0 = 'IGAP'
  ancestry1 = c('UKBB')
  other = paste(ancestry1, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/AD_different/annot_pvalus/result/", ancestry0, "_result_final_10_median_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
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
      # lambdas1 <- append(lambdas1, list(result$lambda_s[[i]]))
    }
    lambdas1 <- append(lambdas1, list(result$lambda_s[[i]]))
  }
  
  
  ancestry0 = 'UKBB'
  ancestry1 = c('IGAP')
  other = paste(ancestry1, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/AD_different/annot_pvalus/result/", ancestry0, "_result_final_10_median_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
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
      # lambdas2 <- append(lambdas2, list(result$lambda_s[[i]]))
    }
    lambdas2 <- append(lambdas2, list(result$lambda_s[[i]]))
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





library(ggplot2)
df <- data.frame(lambdas = as.numeric(lambdas1))
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



df <- data.frame(lambdas = as.numeric(lambdas2))
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,  
    fill = "#0072B2", 
    color = "white",
    alpha = 0.8     
  ) +
  labs(
    title = "Distribution of lambda_s of EUR p-values for AFR GWAS",
    x = "Lambda_s",            
    y = "Frequency (Count)"         
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot


