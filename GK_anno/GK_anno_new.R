
suppressPackageStartupMessages({
  library(glmnet)
  library(xtune)
  library(SNPknock)
  library(adaptMT)
  library(splines)
  library(knockoff)
  library(dplyr)    
  library(corpcor)
  library(MASS)
  library(gam)
  library(randomForest)
  library(tidyverse) 
  library(mgcv)
  library(expm)
  library(rdist)
  library(hdi)
  library(ghostbasil)
})


obj_fun_gk <- function(lam, lambda_0, beta, R, n, d = 20, M = 5){
  p <- length(beta) / (M + 1)
  eta <- as.vector(R %*% lam) / d
  max_eta <- max(eta)
  w <- exp(eta - max_eta) * p / sum(exp(eta - max_eta))
  beta_mat      <- matrix(beta, nrow = p, ncol = M + 1)
  combined_beta <- rowSums(abs(beta_mat))
  
  part1 = - (M+1) * sum(log(w)) / n
  part2 <- sum(w * lambda_0 * combined_beta)
  part3 <- 0.5 * sum(lam^2) / n
  return(part1 + part2 + part3)
}


weight_standardized <- function(lam, R, d = 20){
  p = nrow(R)
  eta <- as.vector(R %*% lam) / d
  max_eta = max(eta)
  w <- exp(eta - max_eta) * p / (sum(exp(eta - max_eta)))
  return(w)
}





suppressPackageStartupMessages({
  library(knockoff)
  library(dplyr)    
})


sol_beta_and_lam_ps = function(fit.prelim, Zscore_0, M, R, temp.A, r_all, n,
                               lambda.seq, lambda_0, init_scale = NULL, maxiter = 100, 
                               verbose = TRUE, init_lam = NULL) {
  # calculate beta if lambda_0 is given, and then update lambda_0 
  if(is.null(init_scale)){
    init_scale = 0
  }
  
  p = length(Zscore_0)
  Prior = numeric((M+1)*p)+1
  
  fit_basil <- ghostbasil(
    temp.A, r_all,
    user.lambdas        = lambda.seq* init_scale,
    alpha               = 1,
    penalty= Prior,
    delta.strong.size   = max(1, min(500, length(r_all)/20)),
    max.strong.size     = length(r_all),
    n.threads           = 1,
    use.strong.rule     = FALSE
  )
  
  beta <- fit_basil$betas[, ncol(fit_basil$betas)]  
  
  
  beta_old = beta
  lambda_s_old = numeric(ncol(R)) + 100
  
  
  lambda_s_all <- numeric()
  objs = c()
  
  init_lam <- numeric(ncol(R)) + 5
  
  for(t in 1:100){
    
    res <- optim(
      par   = init_lam,
      fn    = obj_fun_gk,   
      method = "BFGS", 
      lambda_0 = lambda_0,
      beta  = beta,
      R     = R,
      n     = n,
      M = M
    )
    
    lambda_s = res$par
    print(lambda_s)
    lambda_s_all <- append(lambda_s_all, lambda_s)
    
    w = weight_standardized(lambda_s, R)
    
    Prior = rep(w,M+1)
    
    fit_basil <- ghostbasil(
      temp.A, r_all,
      user.lambdas        = lambda.seq* sum(Prior) / ((M+1)*p),
      alpha               = 1,
      penalty= Prior,
      delta.strong.size   = max(1, min(500, length(r_all)/20)),
      max.strong.size     = length(r_all),
      n.threads           = 1,
      use.strong.rule     = FALSE
    )
    
    beta      <- fit_basil$betas[, ncol(fit_basil$betas)]
    
    
    obj_curr = obj_fun_gk(lambda_s, lambda_0, beta, R, n, M = M)
    obj_old = obj_fun_gk(lambda_s_old, lambda_0, beta_old, R, n, M = M)
    obj_value = obj_curr - obj_old
    objs = c(objs, obj_curr)
    
    
    cat("iter = ", t, " lambda_s = ", lambda_s, " diff lambda_s = ", abs(lambda_s - lambda_s_old), " obj diff = ", obj_value, " sum(Prior) = ", sum(Prior), "\n")
    
    if(mean(abs(lambda_s - lambda_s_old)) < 0.001 || abs(obj_value) < 1e-4){
      print(paste0("convergence achieved at ", t, "-th iteration"))
      break
    }
    
    lambda_s_old = lambda_s
    
    beta_old = beta
  }
  
  list(beta = beta, objs = objs, lambda_s_all = lambda_s_all, Prior = Prior, lambda_s = lambda_s)
}




cv_func <- function(r_all, temp.A, A.left, lambda.seq, Prior, n, nfold = 5){
  nA<-n*(nfold-1)/nfold;nB<-n/nfold
  temp.left<-sqrt(nB/nA/n)*as.matrix(A.left)
  r_all_A<-r_all+as.matrix(temp.left%*%matrix(rnorm(ncol(temp.left)),ncol(temp.left),1))
  r_all_B<-(r_all*n-r_all_A*nA)/nB
  
  fit.basil<-ghostbasil(temp.A, r_all_A, user.lambdas=lambda.seq, alpha=1, delta.strong.size = max(1,min(500,length(r_all_A)/20)), 
                        penalty= Prior, max.strong.size = nrow(temp.A),n.threads=1,use.strong.rule=F)
  beta<-fit.basil$betas[,ncol(fit.basil$betas)]
  # print(hist(beta))
  if(all(beta == 0) == TRUE){
    return(-Inf)
  }
  Get.f<-function(x){x<-as.matrix(x);return(t(x)%*%r_all_B/sqrt(t(x)%*%temp.A%*%x))}
  return(Get.f(beta))
}




GK_anno <- function(Z, R, M = 1, LD, n, ts = 'lasso'){
  
  p = nrow(LD)
  var_para <- apply(R, 2, var)
  
  if(all(is.na(R)) | sum(var_para!=0) == 0){
    print("no informative annotations")
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD,
      M       = M,
      method  = "sdp" 
    )
    GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method='lasso')
    # GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
    T_0 = GK1_lasso$T_0[[1]]
    T_k = GK1_lasso$T_k[[1]]
    
    return(list(T_0 = T_0, T_k = T_k, CV_BEST = NA, lambda_BEST = NA,
                obj_cor_list = NA, lambda_s = rep(0, ncol(R)), Prior = rep(0, 2*p)))
  }
  
  R <- R[, var_para!=0]
  R = scale(R)
  
  # set.seed(123)
  
  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp" 
  )
  
  GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method=ts)
  temp.A = GK1_lasso$temp.A
  r_all = GK1_lasso$r_all
  lambda.seq = GK1_lasso$lambda.seq
  lambda = GK1_lasso$lambda
  
  coeff_list <- exp(seq(-2, 2, by = 0.2))
  
  Prior_final = numeric(2 * p) + 1
  CV_BEST = -Inf
  lambda_BEST = 0
  lambda_s_final = NA
  beta_final = NA
  obj_cor_list <- numeric()
  
  for(coeff in coeff_list){
    
    GKanno_M1_anno_ps <- sol_beta_and_lam_ps(fit.prelim = fit.prelim, 
                                             Zscore_0 = Z, 
                                             M = M, 
                                             R = R, 
                                             temp.A = temp.A, 
                                             r_all = r_all, 
                                             n = n,
                                             lambda.seq = coeff*lambda.seq, 
                                             lambda_0 = coeff*lambda)
    
    obj_cor <- cv_func(r_all = r_all, temp.A = temp.A, A.left = fit.prelim$A.left, lambda.seq = coeff*lambda.seq, 
                       Prior = GKanno_M1_anno_ps$Prior, n = n, nfold = 5)
    
    obj_cor_list <- append(obj_cor_list, obj_cor)
    
    if(obj_cor > CV_BEST){
      lambda_BEST = coeff*lambda
      CV_BEST = obj_cor
      Prior_final = GKanno_M1_anno_ps$Prior
      lambda_s_final = GKanno_M1_anno_ps$lambda_s
      beta_final = GKanno_M1_anno_ps$beta
    }
  }
  
  T_0<-abs(beta_final[1:p])
  T_k<-abs(matrix(beta_final[-(1:p)],p,M))
  
  return(list(T_0 = T_0, T_k = T_k, CV_BEST = CV_BEST, lambda_BEST = lambda_BEST,
              obj_cor_list = obj_cor_list, lambda_s = lambda_s_final, Prior = Prior_final))
}



################# M = 5

GK_anno_M <- function(Z, R, M = 1, LD, n, ts = 'lasso'){
  
  p = nrow(LD)
  var_para <- apply(R, 2, var)
  
  if(all(is.na(R)) | sum(var_para!=0) == 0){
    print("no informative annotations")
    fit.prelim <- GhostKnockoff.prelim(
      cor.G   = LD,
      M       = M,
      method  = "sdp" 
    )
    GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method='lasso')
    # GK.filter<-GhostKnockoff.filter(GK1_lasso$T_0[[1]],GK1_lasso$T_k[[1]])
    T_0 = GK1_lasso$T_0[[1]]
    T_k = GK1_lasso$T_k[[1]]
    
    return(list(T_0 = T_0, T_k = T_k, CV_BEST = NA, lambda_BEST = NA,
                obj_cor_list = NA, lambda_s = rep(0, ncol(R)), Prior = rep(0, 2*p)))
  }
  
  R <- R[, var_para!=0]
  R = scale(R)
  
  # set.seed(123)
  
  fit.prelim <- GhostKnockoff.prelim(
    cor.G   = LD,
    M       = M,
    method  = "sdp" 
  )
  
  GK1_lasso <- GhostKnockoff.fit(Z, n, fit.prelim, method=ts)
  temp.A = GK1_lasso$temp.A
  r_all = GK1_lasso$r_all
  lambda.seq = GK1_lasso$lambda.seq
  lambda = GK1_lasso$lambda
  
  coeff_list <- exp(seq(-2, 2, by = 0.2))
  
  Prior_final = numeric((M+1) * p) + 1
  CV_BEST = -Inf
  lambda_BEST = 0
  lambda_s_final = NA
  beta_final = NA
  obj_cor_list <- numeric()
  
  for(coeff in coeff_list){
    
    GKanno_M1_anno_ps <- sol_beta_and_lam_ps(fit.prelim = fit.prelim, 
                                             Zscore_0 = Z, 
                                             M = M, 
                                             R = R, 
                                             temp.A = temp.A, 
                                             r_all = r_all, 
                                             n = n,
                                             lambda.seq = coeff*lambda.seq, 
                                             lambda_0 = coeff*lambda)
    
    obj_cor <- cv_func(r_all = r_all, temp.A = temp.A, A.left = fit.prelim$A.left, lambda.seq = coeff*lambda.seq, 
                       Prior = GKanno_M1_anno_ps$Prior, n = n, nfold = 5)
    
    obj_cor_list <- append(obj_cor_list, obj_cor)
    
    if(obj_cor > CV_BEST){
      lambda_BEST = coeff*lambda
      CV_BEST = obj_cor
      Prior_final = GKanno_M1_anno_ps$Prior
      lambda_s_final = GKanno_M1_anno_ps$lambda_s
      beta_final = GKanno_M1_anno_ps$beta
    }
  }
  
  T_0<-abs(beta_final[1:p])
  T_k<-abs(matrix(beta_final[-(1:p)],p,M))
  
  return(list(T_0 = T_0, T_k = T_k, CV_BEST = CV_BEST, lambda_BEST = lambda_BEST,
              obj_cor_list = obj_cor_list, lambda_s = lambda_s_final, Prior = Prior_final))
}

