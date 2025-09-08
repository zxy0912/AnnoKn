source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")
packageVersion("Matrix") 




library(tidyr)
library(dplyr)
library(genio)

bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


maf_cal <- function(x){
  freq <- sum(x, na.rm = TRUE)/(2 * sum(!is.na(x)))
  maf <- ifelse(freq < 0.5, freq, 1 - freq)
  return(maf)
}



ancestry = 'EUR'
chrid = '1'
path <- paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/",ancestry,"_0.05maf/bychr/chr_",chrid,"_2w_new")
genotype <- read_plink(path)


# select the X from a LD block:
ld_block <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/EUR/fourier_ls-all.bed", header = TRUE)
ld_block = ld_block[ld_block$chr == 'chr1',]
i = 1

#################### find proper regions

snps <- numeric()
snp_clusters <- numeric()
for(i in 1: nrow(ld_block)){
  
  print(i)
  
  start = ld_block[i,2]
  end = ld_block[i,3]
  index = which(genotype$bim$pos >= start & genotype$bim$pos <= end)
  
  X = t(genotype$X[index, 1:n])
  maf <- apply(X, 2, maf_cal)
  X = bedNA(X)
  
  snps <- append(snps, ncol(X))
  
  #### use hierarchical clustering to cluster SNPs
  
  R = cor(X[1:1000,])
  d <- as.dist(1 - R)  # higher distance = lower correlation
  hc <- hclust(d, method = "average")
  plot(hc, labels = FALSE, main = "SNP clustering by LD")
  clusters <- cutree(hc, h = 0.4)  # e.g., cluster SNPs with r² > 0.6
  length(table(clusters))
  
  snp_clusters = append(snp_clusters, length(table(clusters)))
}

snp_number <- data.frame(index = 1:nrow(ld_block),
                         snps = snps,
                         snp_clusters = snp_clusters)
write.table(snp_number, 
            "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/EUR_0.05maf/real_simulation/region_information/chr1.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

good_iter <- which(snp_number$snp_clusters > 200)





######### start the iteration simulations:

power1_all <- power2_all <- power8_all <- power10_all <- numeric()
fdr1_all <- fdr2_all <- fdr8_all <- fdr10_all <- numeric()

for(iter in 1:50){
  
  i = sample(good_iter, 1)
  
  
  start = ld_block[i,2]
  end = ld_block[i,3]
  index = which(genotype$bim$pos >= start & genotype$bim$pos <= end)
  
  X = t(genotype$X[index, 1:n])
  maf <- apply(X, 2, maf_cal)
  X = bedNA(X)
  
  #### use hierarchical clustering to cluster SNPs
  
  R = cor(X[1:1000,])
  d <- as.dist(1 - R)  # higher distance = lower correlation
  hc <- hclust(d, method = "average")
  plot(hc, labels = FALSE, main = "SNP clustering by LD")
  clusters <- cutree(hc, h = 0.4)  # e.g., cluster SNPs with r² > 0.6
  length(table(clusters))
  mean(table(clusters)) # average size for each cluster
  
  
  ######### select representative SNP
  
  df <- data.frame(
    SNP = names(clusters),
    maf = maf,
    cluster = clusters,
    index = 1:length(maf)
  )
  
  representatives <- df %>%
    group_by(cluster) %>%
    slice_max(order_by = maf, n = 1, with_ties = FALSE)
  
  head(representatives)
  dim(representatives)
  index <- representatives$index[order(representatives$index)]
  X = X[, index]
  
  
  ##########
  
  
  
  n = 10000 # sample size
  # p = 200 # number of variable
  h2e = 0.5 # heritability
  alphae = 0.1 # proportion of causal SNPs
  # rho = 0.5 # maximum correlation between SNPs
  alphalist  <- seq(0.4, 0.05, by = -0.05)
  fdr_seq  <- seq(0.4, 0.05, by = -0.05)
  len = length(fdr_seq)
  

  
  
  X = scale(X)
  p = nrow(representatives)
  betae <- rep(0,p)
  rand <- sample(p,alphae*p)
  
  sigprob = rep(0,p)
  sigprob[1:50] = 1/(1:50)^2/(sum(1/(1:50)^2))
  rand = sample(1:p,alphae*p,prob = sigprob)
  
  print(rand)
  betae[rand] <- rnorm(length(rand),0,sqrt(h2e/p/alphae))
  y <- X%*%betae+rnorm(n,0,sqrt(1-h2e))
  
  
  ss_b = ss_se = ss_z = numeric(p)
  
  for(l in 1:p){
    linearmodel <- lm(y ~ X[,l])
    ss_b[l] <- coef(summary(linearmodel))[2,1]
    ss_se[l] <- coef(summary(linearmodel))[2,2]
    ss_z[l] <- coef(summary(linearmodel))[2,3]
  }
  
  Z       <- ss_z  
  LD <- cor(X)
  n0 = length(rand)
  
  
  ##################
  
  power1 = numeric(len) ## knockoff
  power2 = numeric(len) ## Ghostknockoff M=1 pseudo sum
  power8 = numeric(len) ## GK-anno pseudo sum
  power10 = numeric(len)
  
  
  fdr1 = numeric(len) ## knockoff
  fdr2 = numeric(len) ## Ghostknockoff M=1 pseudo sum
  fdr8 = numeric(len) ## GK-anno pseudo sum
  fdr10 = numeric(len)
  
  
  ############### 1. knockoff ############
  start <- Sys.time()
  
  X_tilde = create.gaussian(X,rep(0,p), LD) # generate knockoff variable
  X_comb = cbind(X,X_tilde)
  
  mdl = cv.glmnet(X_comb,y,alpha=1)
  cvlambda = mdl$lambda.min
  beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
  T0 = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T0,T_tilde)
  W1 = T0-T_tilde
  
  for (j in 1:len) {
    alpha = alphalist[j]
    tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
    rej1 = as.numeric(which(W1>=tau))
    power1[j] = power_cal(rej1, rand)
    fdr1[j] = fdr_cal(rej1, rand)
    
  }
  
  end <- Sys.time()
  cat("knockoff time:", as.numeric(end - start), "seconds\n")
  
  ############# 2. GK M=1 pseudo sum
  start <- Sys.time()
  LD = cor(X)
  N.effect = n
  
  M <- 1
  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp" 
  )
  GK1_lasso <- GhostKnockoff.fit(Z ,N.effect,fit.prelim,method='lasso')
  GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
  
  for (j in 1:len){
    threshold = alphalist[j]
    result_gk <- which(GK.filter$q <= threshold)
    power_gk = length(intersect(result_gk,rand))/n0
    fdp_gk = length(setdiff(result_gk,rand))/max(length(result_gk),1)
    
    print(power_gk)
    print(fdp_gk)
    
    power2[j] <- power_gk
    fdr2[j] <- fdp_gk
  }
  
  end <- Sys.time()
  cat("GK M=1 pseudo sum time:", as.numeric(end - start), "seconds\n")

  
  ############### 8. Ghostknockoff M=1 pseudo sum with annotation
  
  Z       <- ss_z
  z <-1:p
  R <- scale(as.matrix(z))
  M <- 1
  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp" 
  )
  
  GK1_M1_anno_ps <- sol_beta_and_lam_simple_ps(fit.prelim, Z, M, R, GK1_lasso$temp.A, GK1_lasso$r_all, GK1_lasso$lambda.seq, GK1_lasso$lambda, init_scale = NULL, maxiter = 100, 
                                               verbose = TRUE, init_lam = NULL)
  n.G<-length(Z)
  beta <- GK1_M1_anno_ps$beta
  T_0<-abs(beta[1:n.G])
  T_k<-abs(matrix(beta[-(1:n.G)],n.G,M))
  
  GK.filter<-GhostKnockoff.filter(T_0,T_k)
  
  for (j in 1:len){
    threshold = alphalist[j]
    result_gk <- which(GK.filter$q <= threshold)
    power_gk = length(intersect(result_gk,rand))/n0
    fdp_gk = length(setdiff(result_gk,rand))/max(length(result_gk),1)
    
    print(power_gk)
    print(fdp_gk)
    
    power8[j] <- power_gk
    fdr8[j] <- fdp_gk
  }
  

  
  ############### 10. GK-anno M=1 pseudo sum 
  
  Z       <- ss_z
  z <-1:p
  R <- scale(as.matrix(z))
  M <- 1
  
  result = GK_anno(Z, R, M, LD, N.effect)
  
  
  n.G<-length(Z)
  beta <- result$beta_final
  # beta = beta_final
  T_0<-abs(beta[1:n.G])
  T_k<-abs(matrix(beta[-(1:n.G)],n.G,M))
  
  GK.filter<-GhostKnockoff.filter(T_0,T_k)
  
  for (j in 1:len){
    threshold = alphalist[j]
    result_gk <- which(GK.filter$q <= threshold)
    power_gk = length(intersect(result_gk,rand))/n0
    fdp_gk = length(setdiff(result_gk,rand))/max(length(result_gk),1)
    
    print(power_gk)
    print(fdp_gk)
    
    power10[j] <- power_gk
    fdr10[j] <- fdp_gk
  }
  
  print(power2)
  print(power8)
  print(power10)
  
  power1_all <- rbind(power1_all, power1)
  fdr1_all <- rbind(fdr1_all, fdr1)
  
  power2_all <- rbind(power2_all, power2)
  fdr2_all <- rbind(fdr2_all, fdr2)
  
  power8_all <- rbind(power8_all, power8)
  fdr8_all <- rbind(fdr8_all, fdr8)
  
  power10_all <- rbind(power10_all, power10)
  fdr10_all <- rbind(fdr10_all, fdr10)
  
}

colMeans(power1_all)
colMeans(power2_all)
colMeans(power8_all)
colMeans(power10_all)


colMeans(fdr1_all)
colMeans(fdr2_all)
colMeans(fdr8_all)
colMeans(fdr10_all)
