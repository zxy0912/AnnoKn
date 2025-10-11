
########### make the manhatten plot
library("qqman")

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/", "GWAS1",".txt")
sum_s_1 <- read.table(path, header = TRUE)
sum_s <- sum_s_1[sum_s_1$p<0.001,]
p <- manhattan(sum_s, chr = "chr", bp = "pos", p = "p", snp = "SNP", # main=paste("Manhattan plot for the first AD GWAS"), 
               col = c("red4","red","orange","gold","darkgreen","forestgreen","yellowgreen","darkcyan","darkblue","royalblue1", "blue","dodgerblue","deepskyblue","skyblue1","purple4","darkmagenta","violetred","hotpink","palevioletred","lightpink","chocolate4","lightgray","gray28"), 
               ylim=c(0,150),chrlabs = NULL,suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08),highlight = NULL, 
               logp = TRUE, annotatePval = NULL, annotateTop = TRUE)
p

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/", "GWAS2",".txt")
sum_s_2 <- read.table(path, header = TRUE)
sum_s_2$p[which(sum_s_2$p < 10^(-140))] = 10^(-140)

sum_s <- sum_s_2[sum_s_2$p<0.001,]
q <- manhattan(sum_s, chr = "chr", bp = "pos", p = "p", snp = "SNP", # main=paste("Manhattan plot for the second AD GWAS"), 
               col = c("red4","red","orange","gold","darkgreen","forestgreen","yellowgreen","darkcyan","darkblue","royalblue1", "blue","dodgerblue","deepskyblue","skyblue1","purple4","darkmagenta","violetred","hotpink","palevioletred","lightpink","chocolate4","lightgray","gray28"), 
               ylim=c(0,150),chrlabs = NULL,suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08),highlight = NULL, 
               logp = TRUE, annotatePval = NULL, annotateTop = TRUE)
q








merged <- merge(sum_s_1[, c("SNP", "p")], 
                sum_s_2[, c("SNP", "p")], 
                by = "SNP", 
                suffixes = c("_1", "_2"))
merged$logp1 <- -log10(merged$p_1)
merged$logp2 <- -log10(merged$p_2)
correlation <- cor(merged$logp1, merged$logp2, use = "complete.obs")


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

region_1 <- numeric()
region_anno_1 <- numeric()
region_2 <- numeric()
region_anno_2 <- numeric()
lambdas1 <- numeric()
lambdas2 <- numeric()

threshold = 0.1
sum1 = 0
sum2 = 0
chrid = '1'
M = 1
seed = '12345'
# seed = '12345' for dss


chrlist = (1:22)
for(chrid in chrlist){
  print(paste0("chrid", chrid))
  study1 = 'GWAS1'
  study2 = c('GWAS2')
  other = paste(study2, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/AD_new/annot_pvalus/result/", study1, "_result_final_10_median_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)
  
  print(length(result$risk_region))
  
  if(length(result$risk_region) == 0){
    next
  }
  
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
    }
    lambdas1 <- append(lambdas1, result$lambda_s[[i]])
  }
}  
 
chrlist = (1:22)
for(chrid in chrlist){
  print(paste0("chrid", chrid)) 
  study1 = 'GWAS2'
  study2 = c('GWAS1')
  other = paste(study2, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/AD_new/annot_pvalus/result/", study1, "_result_final_10_median_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)
  
  print(length(result$risk_region))
  
  if(length(result$risk_region) == 0){
    next
  }
  
  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)
  
  threshold = 0.1
  
  for(i in 1:n_region){
    # print(i)
    a = q_gk[[i]]
    if(sum(a <= threshold) > 0){
      region_2 <- rbind(region_2, as.numeric(append(chrid, result$risk_region[[i]])))
    }else{
      print(FALSE)
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      region_anno_2 <- rbind(region_anno_2, as.numeric(append(chrid, result$risk_region[[i]])))
    }
    lambdas2 <- append(lambdas2, result$lambda_s[[i]])
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
df <- data.frame(lambdas = lambdas1)
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,  
    fill = "#0072B2", 
    color = "white",
    alpha = 0.8     
  ) +
  labs(
    title = "Distribution of lambda_s of GWAS2 p-values for GWAS1",
    x = "Lambda_s",            
    y = "Frequency (Count)"         
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot



df <- data.frame(lambdas = lambdas2)
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,  
    fill = "#0072B2", 
    color = "white",
    alpha = 0.8     
  ) +
  labs(
    title = "Distribution of lambda_s of GWAS1 p-values for GWAS2",
    x = "Lambda_s",            
    y = "Frequency (Count)"         
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot

