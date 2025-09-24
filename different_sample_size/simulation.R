
source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")

iteration = 100
alphalist  <- seq(0.3, 0.05, by = -0.05)
len = length(alphalist)

power4 = matrix(0, nrow = iteration, ncol = len) ## Ghostknockoff M=1 lasso min
power5 = matrix(0, nrow = iteration, ncol = len) ## AnnoGK M = 1
power6 = matrix(0, nrow = iteration, ncol = len) ## AnnoGK DSS M = 1
power7 = matrix(0, nrow = iteration, ncol = len) ## AnnoGK M = 5
power8 = matrix(0, nrow = iteration, ncol = len) ## AnnoGK DSS M = 5



fdr4 = matrix(0, nrow = iteration, ncol = len)
fdr5 = matrix(0, nrow = iteration, ncol = len)
fdr6 = matrix(0, nrow = iteration, ncol = len)
fdr7 = matrix(0, nrow = iteration, ncol = len)
fdr8 = matrix(0, nrow = iteration, ncol = len)



iter = 1
for(iter in 1:iteration){
  
  print(iter)
  seed <- 1000*iter
  
  N = 1000
  p = 300 ## 200 500 1000
  n0 = 30
  amplitude = 3 ## 1 2 4
  rho = 0.5
  
  
  set.seed(seed)
  
  
  
  
  sigprob = rep(0,p)
  sigprob[1:60] = 1/(1:60)^2/(sum(1/(1:60)^2))
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
  
  X = matrix(rnorm(N*p), N ,p)%*%chol(Covariance)
  X = scale(X)
  beta = rep(0,p)
  beta[rand] = amplitude
  beta = beta*rand_sign
  y = X%*%beta + sqrt(N)*rnorm(N) # generate y via linear model 
  ## var(Y) \approx p
  y = (y - mean(y))/sd(y)
  z <-1:p
  R <- scale(as.matrix(z))
  
  
  
  N.effect <- sample(500:N, p)
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
  
  
  
  # n = mean(N.effect)
  # ts = 'lasso'
  # 
  # 
  # fit.prelim <- GhostKnockoff.prelim(
  #   cor.G   = LD,
  #   M       = M,
  #   method  = "sdp" 
  # )
  # 
  # GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method=ts)
  # 
  # # test if the current version of obtained lambda is the same:
  # set.seed(1000)
  # GK1_lasso <- GhostKnockoff.fit(Z, Neff, fit.prelim, method=ts)
  # GK1_lambda <- GhostKnockoff.lambda(Z, Neff, fit.prelim, method=ts)
  # print(GK1_lasso$lambda.seq)
  # print(GK1_lambda$lambda.seq)
  
  
  
  ################## 4. ghostknockff
  
  ss_b = ss_se = ss_z = numeric(p)
  
  for(l in 1:p){
    linearmodel <- lm(y ~ X[,l])
    ss_b[l] <- coef(summary(linearmodel))[2,1]
    ss_se[l] <- coef(summary(linearmodel))[2,2]
    ss_z[l] <- coef(summary(linearmodel))[2,3]
  }
  
  Z       <- ss_z 
  LD = cor(X)
  M = 1
  
  
  set.seed(seed)
  
  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp" 
  )
  GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
  GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej4 <- which(GK.filter$q <= threshold)
    
    power4[iter, j] = power_cal(rej4, rand)
    fdr4[iter, j] = fdr_cal(rej4, rand)
  }
  
  
  ############# 5. GK1 ps with annotation M = 1
  
  set.seed(seed)
  GK1ps_anno = GK_anno_M(Z, R, M = 1, LD, N.median)
  
  
  # beta <- GK1ps_anno$beta_final
  # T_0<-abs(beta[1:p])
  # T_k<-abs(matrix(beta[-(1:p)],p,M))
  T_0 <- GK1ps_anno$T_0
  T_k <- GK1ps_anno$T_k
  
  GK.filter<-GhostKnockoff.filter(T_0,T_k)
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej5 <- which(GK.filter$q <= threshold)
    power5[iter, j] = power_cal(rej5, rand)
    fdr5[iter,j] = fdr_cal(rej5, rand)
  }
  
  
  ############# 6. GK1 ps with annotation M = 5
  
  # set.seed(seed)
  # GK1ps_anno = GK_anno_M(Z, R, M = 5, LD, N.median)
  # 
  # 
  # # beta <- GK1ps_anno$beta_final
  # # T_0<-abs(beta[1:p])
  # # T_k<-abs(matrix(beta[-(1:p)],p,M))
  # T_0 <- GK1ps_anno$T_0
  # T_k <- GK1ps_anno$T_k
  # 
  # GK.filter<-GhostKnockoff.filter(T_0,T_k)
  # 
  # for (j in 1:len){
  #   threshold = alphalist[j]
  #   rej6 <- which(GK.filter$q <= threshold)
  #   power6[iter, j] = power_cal(rej6, rand)
  #   fdr6[iter,j] = fdr_cal(rej6, rand)
  # }
  
  
  ############# 7. GK1 ps with annotation M = 1
  set.seed(seed)
  GK1ps_anno = GK_anno_dss(Z, R, M = 1, LD, N.effect)
  
  
  # beta <- GK1ps_anno$beta_final
  # T_0<-abs(beta[1:p])
  # T_k<-abs(matrix(beta[-(1:p)],p,M))
  T_0 <- GK1ps_anno$T_0
  T_k <- GK1ps_anno$T_k
  
  GK.filter<-GhostKnockoff.filter(T_0,T_k)
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej7 <- which(GK.filter$q <= threshold)
    power7[iter, j] = power_cal(rej7, rand)
    fdr7[iter,j] = fdr_cal(rej7, rand)
  }
  
  
  ############# 8. GK1 ps with annotation M = 5
  
  # set.seed(seed)
  # GK1ps_anno = GK_anno_dss(Z, R, M = 5, LD, N.effect)
  # 
  # 
  # # beta <- GK1ps_anno$beta_final
  # # T_0<-abs(beta[1:p])
  # # T_k<-abs(matrix(beta[-(1:p)],p,M))
  # T_0 <- GK1ps_anno$T_0
  # T_k <- GK1ps_anno$T_k
  # 
  # GK.filter<-GhostKnockoff.filter(T_0,T_k)
  # 
  # for (j in 1:len){
  #   threshold = alphalist[j]
  #   rej8 <- which(GK.filter$q <= threshold)
  #   power8[iter, j] = power_cal(rej8, rand)
  #   fdr8[iter,j] = fdr_cal(rej8, rand)
  # }
  # 
}

