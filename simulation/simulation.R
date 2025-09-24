

source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")

alphalist  <- seq(0.3, 0.05, by = -0.05)
fdr_seq  <- seq(0.3, 0.05, by = -0.05)
len = length(fdr_seq)
iteration = 100

power1 = matrix(0, nrow = iteration, ncol = len) ## knockoff
power2 = matrix(0, nrow = iteration, ncol = len) ## knockoff simple
power3 = matrix(0, nrow = iteration, ncol = len) ## GK M=1 pseudo sum
power4 = matrix(0, nrow = iteration, ncol = len) ## Ghostknockoff M=1 lasso min
power5 = matrix(0, nrow = iteration, ncol = len) ## Ghostknockoff M=5 pseudo sum
power6 = matrix(0, nrow = iteration, ncol = len) ## Ghostknockoff M=5 lasso min
power7 = matrix(0, nrow = iteration, ncol = len) ## GK M=1 pseudo sum with annotation simple
power8 = matrix(0, nrow = iteration, ncol = len) ## GK M=1 lasso min with annotation simple
power9 = matrix(0, nrow = iteration, ncol = len) ## GK M=1 pseudo sum with annotation
power10 = matrix(0, nrow = iteration, ncol = len) ## GK M=1 lasso min with annotation

fdr1 = matrix(0, nrow = iteration, ncol = len)
fdr2 = matrix(0, nrow = iteration, ncol = len)
fdr3 = matrix(0, nrow = iteration, ncol = len)
fdr4 = matrix(0, nrow = iteration, ncol = len)
fdr5 = matrix(0, nrow = iteration, ncol = len)
fdr6 = matrix(0, nrow = iteration, ncol = len)
fdr7 = matrix(0, nrow = iteration, ncol = len)
fdr8 = matrix(0, nrow = iteration, ncol = len)
fdr9 = matrix(0, nrow = iteration, ncol = len)
fdr10 = matrix(0, nrow = iteration, ncol = len)



for(iter in 1:iteration){
  
  print(iter)
  seed <- 100*iter
  
  N.effect = 1000
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
  
  X = matrix(rnorm(N.effect*p), N.effect ,p)%*%chol(Covariance)
  X = scale(X)
  beta = rep(0,p)
  beta[rand] = amplitude
  beta = beta*rand_sign
  y = X%*%beta + sqrt(N.effect)*rnorm(N.effect) # generate y via linear model 
  ## var(Y) \approx p
  y = (y - mean(y))/sd(y)
  
  
  
  ############### 1. knockoff ############
  X_tilde = create_gaussian_s(X,Sigma,s,chol_V) # generate knockoff variable
  X_comb = cbind(X,X_tilde)
  
  start <- Sys.time()
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
    power1[iter, j] = power_cal(rej1, rand)
    fdr1[iter, j] = fdr_cal(rej1, rand)
    
  }
  
  end <- Sys.time()
  cat("knockoff time:", as.numeric(end - start), "seconds\n")
  
  ############### 2. knockoff simple with annotation
  
  z <-1:p
  R <- scale(as.matrix(z))
  
  result = knockoff_simple(X = X, Xk = X_tilde, y = y, R = R)
  beta = result$beta
  
  
  T0 = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T0,T_tilde)
  W2 = T0-T_tilde
  
  for (j in 1:len) {
    alpha = alphalist[j]
    tau = knockoff.threshold(W2,fdr = alpha,offset = 1)
    rej2 = as.numeric(which(W2>=tau))
    power2[iter, j] = power_cal(rej2, rand)
    fdr2[iter, j] = fdr_cal(rej2, rand)
  }
  
  ############### 3. knockoff-anno with annotation
  
  result = knockoff_anno_improved(X = X, Xk = X_tilde, y = y, R = R)
  beta = result$beta
  
  
  T0 = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T0,T_tilde)
  W3 = T0-T_tilde
  
  for (j in 1:len) {
    
    alpha = alphalist[j]
    tau = knockoff.threshold(W3,fdr = alpha,offset = 1)
    rej3 = as.numeric(which(W3>=tau))
    power3[iter, j] = power_cal(rej3, rand)
    fdr3[iter, j] = fdr_cal(rej3, rand)
    
  }
  
  
  
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
  
  
  ############################## 5. GK with marginal M = 1
  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp" 
  )
  GK1_lasso <- GhostKnockoff.fit.original(Z, N.effect, fit.prelim, method='marginal')
  GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej5 <- which(GK.filter$q <= threshold)
    
    power5[iter, j] = power_cal(rej5, rand)
    fdr5[iter, j] = fdr_cal(rej5, rand)
  }
  
  
  ############################## 6. GK with marginal M = 5
  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = 5,
    method  = "sdp" 
  )
  GK1_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
  GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
  beta_lasso <- cbind(GK1_lasso$T_0[[1]], GK1_lasso$T_k[[1]])
  kappa_ko<-apply(abs(beta_lasso),1,which.max)-1 
  tau_ko<-apply(abs(beta_lasso),1,tau_calculation) 
  q_value_ko<-KO_Filter(tau_ko,kappa_ko,M=M)
  
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej6 <- which(GK.filter$q <= threshold)
    
    power6[iter, j] = power_cal(rej6, rand)
    fdr6[iter, j] = fdr_cal(rej6, rand)
  }
  
  
  
  ############### 7. GK_simple pseudo sum 
  
  GK1_M1_anno_ps  = GK_simple(Z = Z, 
                              R = R, 
                              M = 1, 
                              LD = LD,
                              n = N.effect,
                              ts = 'lasso')
  beta <- GK1_M1_anno_ps$beta
  T_0<-abs(beta[1:p])
  T_k<-abs(matrix(beta[-(1:p)],p,M))
  
  GK.filter<-GhostKnockoff.filter(T_0,T_k)
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej7 <- which(GK.filter$q <= threshold)
    power7[iter, j] = power_cal(rej7, rand)
    fdr7[iter, j] = fdr_cal(rej7, rand)
  }
  
  
  ############# 8. GK1 ps with annotation M = 1
  
  GK1ps_anno = GK_anno(Z, R, M, LD, N.effect)
  
  
  # beta <- GK1ps_anno$beta_final
  # T_0<-abs(beta[1:p])
  # T_k<-abs(matrix(beta[-(1:p)],p,M))
  T_0 <- GK1ps_anno$T_0
  T_k <- GK1ps_anno$T_k
  
  GK.filter<-GhostKnockoff.filter(T_0,T_k)
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej8 <- which(GK.filter$q <= threshold)
    power8[iter, j] = power_cal(rej8, rand)
    fdr8[iter,j] = fdr_cal(rej8, rand)
  }
  
  
  ############# 9. GK1 ps with annotation M = 5
  
  GK1ps_anno = GK_anno_M(Z, R, M = 5, LD, N.effect)
  
  
  # beta <- GK1ps_anno$beta_final
  # T_0<-abs(beta[1:p])
  # T_k<-abs(matrix(beta[-(1:p)],p,M))
  T_0 <- GK1ps_anno$T_0
  T_k <- GK1ps_anno$T_k
  
  GK.filter<-GhostKnockoff.filter(T_0,T_k)
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej9 <- which(GK.filter$q <= threshold)
    power9[iter, j] = power_cal(rej9, rand)
    fdr9[iter,j] = fdr_cal(rej9, rand)
  }
  
}

