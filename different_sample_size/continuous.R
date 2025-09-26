
### different sample size

args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
args <- as.numeric(args)

i <- args[1]

source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")



seed <- 1234*i

N = 5000
p = 300
n0 = 30
heri = 0.05
rho = 0.5
alphalist  <- seq(0.4, 0.1, by = -0.05)
fdr_seq  <- seq(0.4, 0.1, by = -0.05)
len = length(fdr_seq)

for (heri in c(0.05, 0.1, 0.2)){
  for (p in c(300, 600, 1000)){
    set.seed(seed)
    power1 = numeric(len) 
    power2 = numeric(len) 
    power3 = numeric(len) 
    power4 = numeric(len) 
    power5 = numeric(len)
    power6 = numeric(len) 
    
    
    fdr1 = numeric(len)
    fdr2 = numeric(len)
    fdr3 = numeric(len)
    fdr4 = numeric(len)
    fdr5 = numeric(len)
    fdr6 = numeric(len)
    
    
    amplitude = sqrt(heri/(1-heri))
    
    sigprob = rep(0,p)
    sigprob[1:60] = 1/(1:60)^2/(sum(1/(1:60)^2))
    #sigprob[1:60] = 1/(1:60)/(sum(1/(1:60)))
    rand = sample(1:p,n0,prob = sigprob)
    #rand = sample(1:p,n0) 
    rand_sign = sample(c(-1,1), size=p, replace=TRUE)
    
    # AR(1) covariates
    Covariance = toeplitz(rho^(0:(p-1)))
    Sigma = Covariance
    Sigma_inv = solve(Sigma)
    s = create.solve_sdp(Sigma)
    Sigma_0 = matrix(0,2*p,2*p)
    Sigma_0[1:p,1:p] = Sigma
    Sigma_0[1:p,(p+1):(2*p)] = Sigma - diag(s)
    Sigma_0[(p+1):(2*p),1:p] = Sigma - diag(s)
    Sigma_0[(p+1):(2*p),(p+1):(2*p)] = Sigma
    Sigma_0_inv = solve(Sigma_0)
    chol_Sigma_0 = chol(Sigma_0)
    D = diag(s)
    P = diag(p) - solve(Sigma,D)
    V = 2*D - D%*%solve(Sigma,D)
    chol_V = chol(V)
    
    X = matrix(rnorm(N*p),N,p)%*%chol(Covariance)
    X = scale(X)
    beta = rep(0,p)
    beta[rand] = amplitude
    beta = beta*rand_sign
    y = X%*%beta + sqrt(n0)*rnorm(N) # generate y via linear model
    print(var(X%*%beta)/var(y))
    h2e = var(X%*%beta)/var(y)
    y = (y - mean(y))/sd(y)
    
    X_tilde = create_gaussian_s(X,Sigma,s,chol_V) # generate knockoff variable
    X_comb = cbind(X,X_tilde)
    
    
    ################## 1. GhostKnockoff
    
    ss_b = ss_se = ss_z = numeric(p)
    
    N.effect <- sample(100:N, p)
    
    
    linear_weights <- ((100:N) - 99)^2
    N.effect <- sample(100:N, p, prob = linear_weights, replace = FALSE)
    
    
    N.median <- median(N.effect)
    
    
    ss_b = ss_se = ss_z = numeric(p)
    
    for(l in 1:p){
      linearmodel <- lm(y[1:N.effect[l]] ~ X[1:N.effect[l],l])
      ss_b[l] <- coef(summary(linearmodel))[2,1]
      ss_se[l] <- coef(summary(linearmodel))[2,2]
      ss_z[l] <- coef(summary(linearmodel))[2,3]
    }
    
    Z       <- ss_z 
    LD = cor(X)
    M = 1
    
    
    start <- Sys.time()
    set.seed(seed)
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD,
      M       = M,
      method  = "sdp" 
    )
    GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
    GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
    
    end <- Sys.time()
    time1 <- as.numeric(difftime(end, start, units = "secs"))
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej1 <- which(GK.filter$q <= threshold)
      
      power1[j] = power_cal(rej1, rand)
      fdr1[j] = fdr_cal(rej1, rand)
    }
    
    
    ############# 2. AnnoGK
    
    z <-1:p
    R <- scale(as.matrix(z))
    
    
    start <- Sys.time()
    set.seed(seed)
    GK1ps_anno = GK_anno_M(Z, R, M, LD, N.median)
    
    
    beta <- GK1ps_anno$beta_final
    T_0<-GK1ps_anno$T_0
    T_k<-GK1ps_anno$T_k
    
    GK.filter<-GhostKnockoff.filter(T_0,T_k)
    
    end <- Sys.time()
    time2 <- as.numeric(difftime(end, start, units = "secs"))
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej2 <- which(GK.filter$q <= threshold)
      power2[j] = power_cal(rej2, rand)
      fdr2[j] = fdr_cal(rej2, rand)
    }
    
    
    
    ############# 3. AnnoGK with different sample sizes:
    
    start <- Sys.time()
    set.seed(seed)
    GK1ps_anno = GK_anno_dss(Z, R, M, LD, N.effect)
    
    T_0_dss<-GK1ps_anno$T_0
    T_k_dss<-GK1ps_anno$T_k
    
    GK.filter<-GhostKnockoff.filter(T_0_dss,T_k_dss)
    
    end <- Sys.time()
    
    time3 <- as.numeric(difftime(end, start, units = "secs"))
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej3 <- which(GK.filter$q <= threshold)
      power3[j] = power_cal(rej3, rand)
      fdr3[j] = fdr_cal(rej3, rand)
    }
    
    
    ############# 4. ghosknockoff M = 5
    start <- Sys.time()
    M = 5
    set.seed(seed)
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD,
      M       = M,
      method  = "sdp" 
    )
    GK5_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
    GK.filter<-GhostKnockoff.filter(GK5_lasso$T_0[[1]],GK5_lasso$T_k[[1]])
    
    end <- Sys.time()
    time4 <- as.numeric(difftime(end, start, units = "secs"))
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej4 <- which(GK.filter$q <= threshold)
      
      power4[j] = power_cal(rej4, rand)
      fdr4[j] = fdr_cal(rej4, rand)
    }
    
    
    ############# 5. AnnoGK with annotation M = 5
    start <- Sys.time()
    set.seed(seed)
    GK5ps_anno = GK_anno_M(Z, R, M, LD, N.median)
    
    
    beta <- GK5ps_anno$beta_final
    T_0<-GK5ps_anno$T_0
    T_k<-GK5ps_anno$T_k
    
    GK.filter<-GhostKnockoff.filter(T_0,T_k)
    
    end <- Sys.time()
    time5 <- as.numeric(difftime(end, start, units = "secs"))
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej5 <- which(GK.filter$q <= threshold)
      power5[j] = power_cal(rej5, rand)
      fdr5[j] = fdr_cal(rej5, rand)
    }
    
    
    ############# 6. AnnoGK-dss with annotation M = 5
    start <- Sys.time()
    set.seed(seed)
    GK5ps_anno = GK_anno_dss(Z, R, M, LD, N.effect)
    
    
    beta <- GK5ps_anno$beta_final
    T_0<-GK5ps_anno$T_0
    T_k<-GK5ps_anno$T_k
    
    GK.filter<-GhostKnockoff.filter(T_0,T_k)
    
    end <- Sys.time()
    time6 <- as.numeric(difftime(end, start, units = "secs"))
    
    for (j in 1:len){
      threshold = alphalist[j]
      rej6 <- which(GK.filter$q <= threshold)
      power6[j] = power_cal(rej6, rand)
      fdr6[j] = fdr_cal(rej6, rand)
    }
    
    
    
    result = list(power1 = power1, 
                  power2 = power2, 
                  power3 = power3,
                  power4 = power4, 
                  power5 = power5, 
                  power6 = power6, 
                  fdr1 = fdr1,
                  fdr2 = fdr2,
                  fdr3 = fdr3,
                  fdr4 = fdr4,
                  fdr5 = fdr5,
                  fdr6 = fdr6,
                  time1 = time1,
                  time2 = time2,
                  time3 = time3,
                  time4 = time4,
                  time5 = time5,
                  time6 = time6,
                  h2e = h2e
    )
    
    # path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/ghostknockoff/simulation/different_sample_size/result/heri_noneven_",heri,"_n_",N,"_p_",p, "_", i,".RData")
    path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/ghostknockoff/simulation/different_sample_size/result/heri_noneven_final_",heri,"_n_",N,"_p_",p, "_", i,".RData")
    # noneven: different sample size
    # final: use part2 + part3, and also non-standardized weights
    
    dir <- dirname(path)
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    
    
    save(result, file = path)
    
    print(path)
    
  }
  
} 

