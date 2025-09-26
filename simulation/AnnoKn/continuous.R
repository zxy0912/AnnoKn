

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
library(hdi)


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
# source("/home/xz527/Rcode/knockoff_side/knockoff_side/issue_1/knockoff_anno.R")


i = 1

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


power1 = numeric(len)
power2 = numeric(len)
power3 = numeric(len)
power4 = numeric(len)
power5 = numeric(len)
power6 = numeric(len)
power7 = numeric(len)
power8 = numeric(len)
fdr1 = numeric(len)
fdr2 = numeric(len)
fdr3 = numeric(len)
fdr4 = numeric(len)
fdr5 = numeric(len)
fdr6 = numeric(len)
fdr7 = numeric(len)
fdr8 = numeric(len)

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


# cvfit <- cv.glmnet(X, y)
# lambda_0 <- cvfit$lambda.min
# Prior = (1:p)
# Prior = Prior*p/sum(Prior)
# beta <- as.vector(glmnet(X, y, lambda = lambda_0, penalty.factor = Prior, intercept = FALSE, standardize = FALSE)$beta)


####################################
## Computing p values
####################################
mdl = lm(y~X)
pvals = summary(mdl)$coefficients[-1,4]

################ 1. AdaPT
z <-1:p
if(binary == 'binary'){
  z = ifelse(z <= 450, 0 , 1) 
}





############# 2. Vanilla  knockoff
Xk = knockoffHMM(X, pInit, Q,pEmit)
mdl = cv.glmnet(cbind(X,Xk),y,alpha=1)
cvlambda = mdl$lambda.min
beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
T = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T,T_tilde)
W1 = T-T_tilde

for (j in 1:len) {
  alpha = alphalist[j]
  tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
  rej2 = as.numeric(which(W1>=tau))
  power2[j] = power_cal(rej2, nonzero)
  fdr2[j] = fdr_cal(rej2, nonzero)
  
}


######### 7. Knockoff-anno
R <- scale(as.matrix(z))
start = Sys.time()
# result_raw <- knockoff_anno_raw(X = X, Xk = Xk, y = y, attempts = c(0), R = R)
result_raw <- knockoff_anno_improved(X = X, Xk = Xk, y = y, attempts = c(0), R = R)
end = Sys.time()
time = as.numeric(difftime(end, start, units = "secs"))

cv_matrix = result_raw$cv_error

lambda_vals <- as.numeric(unlist(result_raw$lambda_s))

for (j in 1:len) {
  
  alpha = alphalist[j]
  acc = calc_power_fdr(result_raw$beta, p, alpha, nonzero)
  power7[j] = acc[1]
  fdr7[j] = acc[2]
  
}

result_raw_old = result_raw






plot(result_raw$lams, result_raw$cv_error, 
     xlab = "Lambda (λ)",       
     ylab = "CV Error",          
     main = "CV Error vs. Lambda Comparison", 
     col = "blue")              
points(result_raw_old$lams, result_raw_old$cv_error, 
       pch = 3, col = "red")              
legend("bottomright",             
       legend = c("Non-standardized", "Standardized"), 
       col = c("blue", "red"),  
       lty = c(1, -1),         
       pch = c(-1, 3))


