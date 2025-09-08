
source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")
packageVersion("Matrix") 

chrid = '3'
ancestry = 'EUR'
maf = 0.01
M = 1

library(tidyr)
library(dplyr)
library(genio)
library(bigstatsr)
library(dendextend)

bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


################ load the summary statistics

# path = paste("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/real_data/ancestry/",ancestry,"/T2D2024/t2d_ss_2024", sep="")
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/T2D/", ancestry, "_Metal_LDSC-CORR_Neff.v2.txt")
print(path)
t2d_ss_2024_0 <- read.table(path, header = T)
t2d_ss_2024 = t2d_ss_2024_0[t2d_ss_2024_0$Chromsome == chrid,]


################ load the reference panel from 1KG

# if(maf == 0.05){
#   maf_indi = 'maf5'
#   maf_dir = 'bychr_maf0.05'
# }else if(maf == 0.01){
#   maf_indi = 'maf1'
#   maf_dir = 'bychr_maf0.01'
# }
# path <- paste("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/1KG_sum/1000G_hg37/EUR/bychr", "/1KG_chr", chrid,sep='')
path <- paste("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/simu_sep/1000G/AFR/bychr_", ancestry, "/1KG_chr", chrid,sep='')
ref_chr <- read_plink(path)

############### load the annotation data

library(data.table)
cats = c("SNP","Coding_UCSC", "TSS_Hoffman", "Promoter_UCSC", "Intron_UCSC", "UTR_3_UCSC", "UTR_5_UCSC")
annot = fread(paste0("/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/baselineLD_2.1/baselineLD.", chrid, ".annot.gz"))
annot = annot[, ..cats]


############### find the common part of three datasets:

common_id <- intersect(ref_chr$bim$id, annot$SNP)
index1 = match(common_id, ref_chr$bim$id)
index2 = match(common_id, annot$SNP)

bim = ref_chr$bim[index1, ]
X = ref_chr$X[index1, ]
table(rownames(X) == bim$id)
annot = annot[index2, ]
table(rownames(X) == annot$SNP)

t2d_ss_2024 <- t2d_ss_2024[order(t2d_ss_2024$Position),]
common_pos <- intersect(t2d_ss_2024$Position, bim$pos)

index1 = match(common_pos, t2d_ss_2024$Position)
index2 = match(common_pos, bim$pos)

t2d_ss_2024 <- t2d_ss_2024[index1, ]
bim = bim[index2, ]
X = X[index2, ]
table(rownames(X) == bim$id)
annot = annot[index2, ]
table(rownames(X) == annot$SNP)
table(t2d_ss_2024$Position == bim$pos)



################### define risk regions

len_half = 250000
risk_id <- which(t2d_ss_2024$Pval < 5 * 10^-8)
pos_id = unlist(t2d_ss_2024$Position[risk_id])
pos_diff <- diff(pos_id)
where_sep = as.numeric(c(0, which(pos_diff > 2 * len_half)))
risk_region = numeric()


if(length(where_sep) > 1){
  for(i in 1:(length(where_sep)-1)){
    risk_region <- rbind(risk_region, c(pos_id[where_sep[i] + 1] - len_half, pos_id[where_sep[i+1]] + len_half))
  }
  risk_region <- rbind(risk_region, c(pos_id[where_sep[length(where_sep)] + 1] - len_half, pos_id[length(pos_id)] + len_half))
}else if(length(where_sep) == 1){
  risk_region <- rbind(risk_region, c(pos_id[where_sep[length(where_sep)] + 1] - len_half, pos_id[length(pos_id)] + len_half))
}




################### address allele flipping

q_gk <- list()
q_gkanno <- list()
lambda_s <- list()


for(i in 1:nrow(risk_region)){
  
  
  set.seed(1000)
  
  start = risk_region[i, 1]
  end = risk_region[i, 2]
  
  
  region_id = which(bim$pos > start & bim$pos < end)
  # table(rownames(X) == annot$SNP)
  # table(bim$id == annot$SNP)
  # table(t2d_ss_2024$Position == bim$pos)
  t2d_region <- t2d_ss_2024[region_id, ]
  bim_region <- bim[region_id, ]
  X_region <- bedNA(t(X[region_id, ]))
  annot_region <- annot[region_id, ]
  # table(colnames(X_region) == annot_region$SNP)
  # table(bim_region$id == annot_region$SNP)
  # table(t2d_region$Position == bim_region$pos)
  
  
  z_scores = t2d_region$Beta/t2d_region$SE
  R = cor(X_region)
  nsample <- mean(t2d_region$Neff)
  
  
  library(susieR)
  lambda = estimate_s_rss(z_scores, R, n = nsample)
  print(lambda)
  
  
  inv <- which(t2d_region$EffectAllele==bim_region$ref & t2d_region$NonEffectAllele==bim_region$alt)
  print(length(inv))
  
  t2d_region$Beta[inv] <- -t2d_region$Beta[inv]
  t2d_region$EAF[inv] <- 1-t2d_region$EAF[inv]
  a <- t2d_region$EffectAllele[inv]
  t2d_region$EffectAllele[inv] <- t2d_region$NonEffectAllele[inv]
  t2d_region$NonEffectAllele[inv] <- a
  
  print(table(t2d_region$EffectAllele==bim_region$alt & t2d_region$NonEffectAllele==bim_region$ref))
  
  
  falseid <- which((t2d_region$EffectAllele==bim_region$alt & t2d_region$NonEffectAllele==bim_region$ref)==F)
  
  print(paste("number of falseid:",length(falseid)))
  
  if(length(falseid)>0){
    print("removing false SNPs:")
    print(falseid)
    t2d_region <- t2d_region[-falseid,]
    bim_region <- bim_region[-falseid,]
    X_region <- X_region[-falseid,]
  }
  
  print(table(t2d_region$EffectAllele==bim_region$alt & t2d_region$NonEffectAllele==bim_region$ref))
  
  
  z_scores = t2d_region$Beta/t2d_region$SE
  
  
  # library(susieR)
  lambda = estimate_s_rss(z_scores, R, n = nsample)
  print(lambda)
  
  
  
  ################### perform variable selection
  
  
  #### use hierarchical clustering to cluster SNPs
  
  d <- as.dist(1 - R)  # higher distance = lower correlation
  hc <- hclust(d, method = "average")
  hc <- as.dendrogram(hc)
  clusters <- cutree(hc, h = 0.4)  # e.g., cluster SNPs with r² > 0.6
  length(table(clusters))
  mean(table(clusters)) # average size for each cluster
  
  
  ######### select representative SNP
  
  table(names(clusters) == bim_region$id)
  
  df <- data.frame(
    SNP = names(clusters),
    pvalues = t2d_region$Pval,
    cluster = clusters,
    index = 1:nrow(t2d_region)
  )
  
  representatives <- df %>%
    group_by(cluster) %>%
    slice_min(order_by = pvalues, n = 1, with_ties = FALSE)
  
  head(representatives)
  dim(representatives)
  index <- representatives$index[order(representatives$index)]
  X_final = X_region[, index]
  z_scores_final  = z_scores[index]
  annot_final = annot_region[index,]
  
  ##########
  X_final = scale(X_final)
  table(colnames(X_final) == annot_final$SNP)
  LD = cor(X_final)
  p = nrow(LD)
  
  
  ############################# ghostknockoff
  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp" 
  )
  GK1_lasso <- GhostKnockoff.fit(z_scores_final, nsample, fit.prelim, method='lasso')
  GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
  threshold = 0.1
  which(GK.filter$q <= threshold)
  a = GK.filter$q
  
  q_gk <- append(q_gk, list(a))
  
  
  ############################# GK-anno
  
  annot_final = annot_final[, 2:ncol(annot_final)]
  annot_final_df <- as.data.frame(annot_final)
  annot_final_df <- annot_final_df[, which(apply(annot_final_df, 2, var) > 0)]
  R_anno = scale(annot_final_df)
  GK1ps_anno = GK_anno(z_scores_final, R_anno, M, LD, nsample)
  
  # beta <- GK1ps_anno$beta_final
  T_0<-GK1ps_anno$T_0
  T_k<-GK1ps_anno$T_k
  
  GK.filter<-GhostKnockoff.filter(T_0,T_k)
  threshold = 0.1
  which(GK.filter$q <= threshold)
  b = GK.filter$q
  
  q_gkanno <- append(q_gkanno, list(b))
  
  lambda_s <- append(lambda_s, list(GK1ps_anno$lambda_s))
  
  
}

threshold = 0.1
sum1 = 0
sum2 = 0
sum3 = 0

for(i in 1:nrow(risk_region)){
  print(i)
  a = q_gk[[i]]
  print(which(a <= threshold))
  if(sum(a <= threshold) > 0){
    sum1 = sum1 + 1
  }
  b = q_gkanno[[i]]
  print(which(b <= threshold))
  if(sum(b <= threshold) > 0){
    sum2 = sum2 + 1
  }
  if(sum(a <= threshold) > 0 & sum(b <= threshold) > 0){
    sum3 = sum3 + 1
  }
}

result <- list(
  q_gk = q_gk,
  q_gkanno = q_gkanno,
  lambda_s = lambda_s)


path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/annot_LDSC/result/result_chr",chrid,".RData")
save(result, file = path)


