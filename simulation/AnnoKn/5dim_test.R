
######################### Simulation for n > p & 5 - dim (1 informative + 4 non-informative) ###########################
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)
args <- as.numeric(args)

i <- args[1]

####################################
library(glmnet)
library("xtune")
library("SNPknock")
suppressMessages(library("adaptMT"))
suppressMessages(library("splines"))
suppressMessages(library("knockoff"))
suppressMessages(library("SNPknock"))
suppressMessages(library("dplyr"))
suppressMessages(library("corpcor"))
suppressMessages(library("glmnet"))
suppressMessages(library("MASS"))
suppressMessages(library("gam"))
suppressMessages(library("randomForest"))
suppressMessages(library("tidyverse"))
suppressMessages(library("mgcv"))
# suppressMessages(library("flare"))
suppressMessages(library("expm"))
suppressMessages(library("rdist"))


source_gitfile <- function(filename){
  source(sprintf("%s.R", filename))
}
file_vec <- c("/home/xz527/Rcode/knockoff_side/code//utils/all_other_methods",
              "/home/xz527/Rcode/knockoff_side/code/utils/adaptive_knockoff",
              "/home/xz527/Rcode/knockoff_side/code/utils/filter_EM",
              "/home/xz527/Rcode/knockoff_side/code/utils/filter_gam",
              "/home/xz527/Rcode/knockoff_side/code/utils/filter_glm",
              "/home/xz527/Rcode/knockoff_side/code/utils/filter_randomForest", 
              "/home/xz527/Rcode/knockoff_side/code/utils/All_q_est_functions", 
              "/home/xz527/Rcode/knockoff_side/code/utils/accumulation_test_functions")
getfile <- sapply(file_vec,source_gitfile)



source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
#source("/home/cl2667/Anno/EM.R")
suppressMessages(library("knockoff"))
suppressMessages(library("SNPknock"))
suppressMessages(library("glmnet"))



# alphalist = seq(0.3,0.01,-0.01)
rho = 0.5

alphalist <- seq(0.3,0.05,-0.05)
len = length(alphalist)


print(paste0("######################################### iteration ", i, "  ######################################### "))

set.seed(100 * i)
n = 1000
p = 900
k = 150

alphalist <- seq(0.3,0.05,-0.05)
Sigma = diag(rep(1,p))
amp = 3.5
binary = ''

for(binary in c('')){
  #for(binary in c('', 'binary')){
  power1 = numeric(len)
  power2 = numeric(len)
  power3 = numeric(len)
  power4 = numeric(len)
  power5 = numeric(len)
  power6 = numeric(len)
  power7 = numeric(len)
  
  fdr1 = numeric(len)
  fdr2 = numeric(len)
  fdr3 = numeric(len)
  fdr4 = numeric(len)
  fdr5 = numeric(len)
  fdr6 = numeric(len)
  fdr7 = numeric(len)
  
  cv_matrix <- numeric(61)
  cv_matrix_em <- numeric(61)
  
  
  Sigma = toeplitz(rho^(0:(p-1)))
  amp = 3.5
  sigprob = rep(0,p)
  sigprob[1:300] = 1/(1:300)^2/(sum(1/(1:300)^2))
  nonzero = sample(1:p,k,prob = sigprob)
  beta0 = amp * (1:p %in% nonzero)*sign(rnorm(p)) / sqrt(n)
  y.sample = function(X) X%*%beta0+rnorm(n,0,1)
  all_res <- data.frame()
  
  ####################################
  ## HMM parameters
  ####################################
  # Number of possible states for each variable
  K=5;M=3;  
  # Marginal distribution for the first variable
  pInit = rep(1/K,K)
  # Create p-1 transition matrices
  Q = array(stats::runif((p-1)*K*K),c(p-1,K,K))
  for(j in 1:(p-1)) { Q[j,,] = Q[j,,] / rowSums(Q[j,,]) }
  pEmit = array(stats::runif(p*M*K),c(p,M,K))
  for(j in 1:p) { pEmit[j,,] = pEmit[j,,] / rowSums(pEmit[j,,]) }
  
  ####################################
  ## Generating data
  ####################################
  X = sampleHMM(pInit, Q, pEmit, n=n)
  y = y.sample(X)
  Xk = knockoffHMM(X, pInit, Q,pEmit)
  
  ################ Informative Annotation ################
  z <-1:p
  if(binary == 'binary'){
    z = ifelse(z <= 450, 0 , 1) 
  }
  R_inf <- as.matrix(z)
  
  ############### non-informative Annotation ##############
  p         <- nrow(R_inf)
  num_inf   <- ncol(R_inf)     # 1
  num_noise <- 4                # dimension of noise
  total     <- num_inf + num_noise
  
  noise_mat <- matrix(
    sample(1:p, size = p * num_noise, replace = TRUE),
    nrow = p, ncol = num_noise
  )
  
  if (binary == "binary") {
    noise_mat <- ifelse(noise_mat <= 450, 0, 1)
  }
  
  causal_anno_index <- sample(seq_len(total), size = num_inf)
  #causal_anno_index <- 1
  
  R5 <- matrix(NA, nrow = p, ncol = total)
  
  R5[, causal_anno_index] <- R_inf
  
  other_idx <- setdiff(seq_len(total), causal_anno_index)
  R5[, other_idx] <- noise_mat
  
  R <- scale(R5)
  cat("Informative annotation index are:", causal_anno_index, "\n")
  
  R_inf <- scale(as.matrix(z))
  
  
  ############# 
  mdl = cv.glmnet(cbind(X,Xk),y,alpha=1)
  cvlambda = mdl$lambda.min
  beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W1 = T-T_tilde
  
  
  ######### 1. Adaknockoff with EM algorithm
  resj =filter_EM(W1,R,alpha =alphalist,offset=1,df = 2)
  for(j in 1:len){
   rej1  = resj$rejs[[j]]
   power1[j] = power_cal(rej1, nonzero)
   fdr1[j] = fdr_cal(rej1, nonzero)
  }

  ######### 2. Adaknockoff with EM algorithm with filter
  resj =filter_EM(W1,R_inf,alpha =alphalist,offset=1,df = 2)
  for(j in 1:len){
   rej2  = resj$rejs[[j]]
   power2[j] = power_cal(rej2, nonzero)
   fdr2[j] = fdr_cal(rej2, nonzero)
  }
  
  ######### 3. Knockoff-anno
  start = Sys.time()
  # result_raw <- knockoff_anno_raw(X = X, Xk = Xk, y = y, attempts = c(0), R = R)
  result_raw <- knockoff_anno_improved(X = X, Xk = Xk, y = y, attempts = c(0), R = R)
  end = Sys.time()
  time = as.numeric(difftime(end, start, units = "secs"))
  
  cv_matrix = result_raw$cv_error
  
  lambda_vals <- as.numeric(unlist(result_raw$lambda_s))
  lambda_inf <- lambda_vals[causal_anno_index]
  lambda_noninf <- lambda_vals[other_idx]
  
  for (j in 1:len) {
    alpha = alphalist[j]
    acc = calc_power_fdr(result_raw$beta, p, alpha, nonzero)
    power3[j] = acc[1]
    fdr3[j] = acc[2]
    
  }
  
  ############# 4. knockoff_anno with filter
  start = Sys.time()
  # result_raw <- knockoff_anno_raw(X = X, Xk = Xk, y = y, attempts = c(0), R = R)
  result_raw <- knockoff_anno_improved(X = X, Xk = Xk, y = y, attempts = c(0), R = R_inf)
  end = Sys.time()
  time = as.numeric(difftime(end, start, units = "secs"))
  
  cv_matrix = result_raw$cv_error
  
  lambda_vals <- as.numeric(unlist(result_raw$lambda_s))
  
  for (j in 1:len) {
    alpha = alphalist[j]
    acc = calc_power_fdr(result_raw$beta, p, alpha, nonzero)
    power4[j] = acc[1]
    fdr4[j] = acc[2]
    
  }
  
  ############# 5. AdaKn(RF) ###################
  resj = filter_randomForest(W1,R,alpha =alphalist,offset=1)
  for(j in 1:len){
    rej5  = resj$rejs[[j]]
    power5[j] = power_cal(rej5, nonzero)
    fdr5[j] = fdr_cal(rej5, nonzero)
  }
  
  
  ############ 6. AdaKn(RF) with filter
  resj = filter_randomForest(W1,R_inf,alpha =alphalist,offset=1)
  for(j in 1:len){
    rej6  = resj$rejs[[j]]
    power6[j] = power_cal(rej6, nonzero)
    fdr6[j] = fdr_cal(rej6, nonzero)
  }
  
  ############### 7. knockoff
  Xk = knockoffHMM(X, pInit, Q,pEmit)
  mdl = cv.glmnet(cbind(X,Xk),y,alpha=1)
  cvlambda = mdl$lambda.min
  beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W7 = T-T_tilde
  
  for (j in 1:len) {
    alpha = alphalist[j]
    tau = knockoff.threshold(W7,fdr = alpha,offset = 1)
    rej7 = as.numeric(which(W7>=tau))
    power7[j] = power_cal(rej7, nonzero)
    fdr7[j] = fdr_cal(rej7, nonzero)
    
  }
  
  
  result = list(power1 = power1, # AdaKn(EM)
                power2 = power2, # AdaKn(EM) with filter
                power3 = power3, # Knockoff Anno
                power4 = power4, # Knockoff Anno with filter
                power5 = power5, # AdaKn(RF)
                power6 = power6, # AdaKn(RF) with filter\
                power7 = power7,
                lambda_inf = lambda_inf,
                lambda_noninf = lambda_noninf,
                fdr1 = fdr1,
                fdr2 = fdr2,
                fdr3 = fdr3,
                fdr4 = fdr4,
                fdr5 = fdr5,
                fdr6 = fdr6,
                fdr7 = fdr7,
                cv_matrix = cv_matrix)
  
  
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/knockoff/simulation/5dimen/result/result_test_", binary, "_amp_",amp, "_", i,".RData")
  
  dir <- dirname(path)
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  if (file.exists(path) && file.info(path)$isdir) {
    unlink(path, recursive = TRUE)
  }
  
  if (file.exists(path) && !file.info(path)$isdir) {
    file.remove(path)
  }
  
  save(result, file = path)
  
  print(path)
  
  
}