








ancestry0 = 'AFR'
sum1_all <- numeric()
sum2_all <- numeric()
idex1 <- numeric()
idex2 <- numeric()

for(chrid in as.character(1)){
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/T2D/annot_pvalus/result/result_", ancestry0, "_chr_EUR",chrid,".RData")
  load(path)
  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)
  
  threshold = 0.1
  sum1 = 0
  sum2 = 0
  
  for(i in 1:n_region){
    a = q_gk[[i]]
    # print(which(a <= threshold))
    if(sum(a <= threshold) > 0){
      sum1 = sum1 + 1
      idex1 = append(idex1, i)
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      sum2 = sum2 + 1
      idex2 = append(idex2, i)
    }
  }
  print(sum1)
  sum1_all <- append(sum1_all, sum1)
  print(sum2)
  sum2_all <- append(sum2_all, sum2)
}


cbind(sum1_all, sum2_all)
sum(sum1_all)
sum(sum2_all)
print(idex1)
print(idex2)









ancestry0 = 'EUR'
sum1_all <- numeric()
sum2_all <- numeric()

for(chrid in as.character(1:22)){
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/T2D/annot_pvalus/result/result_", ancestry0, "_chr_AFR",chrid,".RData")
  load(path)
  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)
  
  threshold = 0.1
  sum1 = 0
  sum2 = 0
  
  for(i in 1:n_region){
    print(i)
    a = q_gk[[i]]
    # print(which(a <= threshold))
    if(sum(a <= threshold) > 0){
      sum1 = sum1 + 1
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      sum2 = sum2 + 1
    }
  }
  print(sum1)
  sum1_all <- append(sum1_all, sum1)
  print(sum2)
  sum2_all <- append(sum2_all, sum2)
}


cbind(sum1_all, sum2_all)
sum(sum1_all)
sum(sum2_all)







chrid = 1

ancestry0 = 'EUR'
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/T2D/annot_pvalus/risk_region/",ancestry0,"_chr_", chrid,".txt")
risk_region1 <- read.table(path)



ancestry0 = 'AFR'
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/T2D/annot_pvalus/risk_region/",ancestry0,"_chr_", chrid,".txt")
risk_region2 <- read.table(path)


common_rows <- merge(risk_region1, risk_region2)
nrow(common_rows)



########## check the risk region for different CHRs
s = numeric()

for(chrid in 1:22){
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/T2D/annot_pvalus/risk_region/","EUR","_chr_", chrid,".txt")
  risk_region <- read.table(path)
  print(dim(risk_region))
  s = append(s, nrow(risk_region))
}



##############


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
M = 3


for(chrid in 1:22){
  print(paste0("chrid", chrid))
  ancestry0 = 'EUR'
  ancestry1 = c('AFR')
  other = paste(ancestry1, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/T2D/annot_pvalus/result/", ancestry0, "_result_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
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
      lambdas1 <- append(lambdas1, result$lambda_s[[i]])
    }
  }
  
  
  ancestry0 = 'AFR'
  ancestry1 = c('EUR')
  other = paste(ancestry1, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/T2D/annot_pvalus/result/", ancestry0, "_result_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)
  print(length(result$risk_region))
  
  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)
  
  threshold = 0.1
  
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
      lambdas2 <- append(lambdas2, result$lambda_s[[i]])
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



nrow(region_1) + nrow(region_2) - 263
nrow(region_anno_1) + nrow(region_anno_2) - 288


