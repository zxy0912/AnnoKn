# library(devtools)
# install_github("GuJQ5/KnockoffCDGE")
library(JuliaCall)
library(KnockoffCDGE)
library(knockoff)
julia <- julia_setup("/home/xz527/julia-1.11.6/bin",installJulia = FALSE)


iteration = 10
alphalist  <- seq(0.3, 0.05, by = -0.05)
len = length(alphalist)
power1 = matrix(0, nrow = iteration, ncol = len) 
power2 = matrix(0, nrow = iteration, ncol = len)
fdr1 = matrix(0, nrow = iteration, ncol = len)
fdr2 = matrix(0, nrow = iteration, ncol = len)


iter = 1
for(iter in 1:iteration){
  
  seed <- 100*iter
  
  N.effect = 1000
  p = 300 ## 200 500 1000
  n0 = 30
  amplitude = 4 ## 1 2 4
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
  
  
  result1 = KnockoffCDGE(Z, Sigma = LD, M = 1, n = N.effect, method = "SDP")
  result5 = KnockoffCDGE(Z, Sigma = LD, M = 5, n = N.effect, method = "SDP")
  
  
  
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej1 <- which(result1$q_value_ko <= threshold)
    
    power1[iter, j] = power_cal(rej1, rand)
    fdr1[iter, j] = fdr_cal(rej1, rand)
  }
  
  for (j in 1:len){
    threshold = alphalist[j]
    rej2 <- which(result5$q_value_ko <= threshold)
    
    power2[iter, j] = power_cal(rej2, rand)
    fdr2[iter, j] = fdr_cal(rej2, rand)
  }
  
}



resultME = rbind(alphalist,
      colMeans(power1),
      colMeans(power2),
      colMeans(fdr1, na.rm = TRUE),
      colMeans(fdr2, na.rm = TRUE))
