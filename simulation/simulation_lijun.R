
# the sample size could be different for different SNPs in summary statistics
# in this simulation, we assume all the sample sizes are the same: N.effect
# however, note that the input for GK_anno_M: N is a constant
# the input for GK_anno_dss and GhostKnockoff: N can be a vector
# based on simulations, GK_anno_dss doesn't perform well compared with GK_anno_M using median as the sample size

source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")



seed <- 1234
set.seed(seed)



N.effect = 5000 # number of samples 
p = 300 # number of SNPs
n0 = 30 # number of causal SNPs
heri = 0.1 # heritability, control the signal strength
rho = 0.5 # AR1 correlation parameter
alphalist  <- seq(0.4, 0.1, by = -0.05)
fdr_seq  <- seq(0.4, 0.1, by = -0.05)
len = length(fdr_seq)




power1 = numeric(len)  # model-X knockoff
power2 = numeric(len)  # AnnoKn-simple
power3 = numeric(len)  # AnnoKn
power4 = numeric(len)  # GhostKnockoff
power5 = numeric(len)  # AnnoGK-simple
power6 = numeric(len)  # AnnoGK
power7 = numeric(len)  # AnnoGK with different sample sizes
power8 = numeric(len)  # Ghosknockoff M = 5
power9 = numeric(len)  # AnnoGK-simple M = 5
power10 = numeric(len) # AnnoGK M = 5
power11 = numeric(len) # AnnoGK with different sample sizes M = 5

fdr1 = numeric(len)
fdr2 = numeric(len)
fdr3 = numeric(len)
fdr4 = numeric(len)
fdr5 = numeric(len)
fdr6 = numeric(len)
fdr7 = numeric(len)
fdr8 = numeric(len)
fdr9 = numeric(len)
fdr10 = numeric(len)
fdr11 = numeric(len)


########### generate the data and use the index as annotation (the probability to be causal is decided by index)
amplitude = sqrt(heri/(1-heri))

sigprob = rep(0,p)
sigprob[1:60] = 1/(1:60)^2/(sum(1/(1:60)^2))
rand = sample(1:p,n0,prob = sigprob)
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

X = matrix(rnorm(N.effect*p),N.effect,p)%*%chol(Covariance)
X = scale(X)
beta = rep(0,p)
beta[rand] = amplitude
beta = beta*rand_sign
y = X%*%beta + sqrt(n0)*rnorm(N.effect) # generate y via linear model
print(var(X%*%beta)/var(y))
h2e = var(X%*%beta)/var(y)
y = (y - mean(y))/sd(y)

X_tilde = create_gaussian_s(X,Sigma,s,chol_V) # generate knockoff variable
X_comb = cbind(X,X_tilde)


############### 1. model-X knockoff ############

start <- Sys.time()
X_tilde = create_gaussian_s(X,Sigma,s,chol_V) # generate knockoff variable
X_comb = cbind(X,X_tilde)

mdl = cv.glmnet(X_comb,y,alpha=1)
cvlambda = mdl$lambda.min
beta = mdl$glmnet.fit$beta[,mdl$lambda ==mdl$lambda.min]
T0 = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T0,T_tilde)
W1 = T0-T_tilde

end <- Sys.time()
time1 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len) {
  alpha = alphalist[j]
  tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
  rej1 = as.numeric(which(W1>=tau))
  power1[j] = power_cal(rej1, rand)
  fdr1[j] = fdr_cal(rej1, rand)
}


############### 2. AnnoKn-simple ###############
z <-1:p
R <- scale(as.matrix(z))

start <- Sys.time()
set.seed(seed)
result = knockoff_simple(X = X, Xk = X_tilde, y = y, R = R)
beta = result$beta
end <- Sys.time()
time2 <- as.numeric(difftime(end, start, units = "secs"))

T0 = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T0,T_tilde)
W2 = T0-T_tilde

for (j in 1:len) {
  alpha = alphalist[j]
  tau = knockoff.threshold(W2,fdr = alpha,offset = 1)
  rej2 = as.numeric(which(W2>=tau))
  power2[j] = power_cal(rej2, rand)
  fdr2[j] = fdr_cal(rej2, rand)
}


############### 3. AnnoKn ###############
start <- Sys.time()
set.seed(seed)
result = knockoff_anno_improved(X = X, Xk = X_tilde, y = y, R = R)
beta = result$beta

end <- Sys.time()
time3 <- as.numeric(difftime(end, start, units = "secs"))


T0 = abs(beta[1:p])
T_tilde = abs(beta[(p+1):(2*p)])
T_max = pmax(T0,T_tilde)
W3 = T0-T_tilde

for (j in 1:len) {
  
  alpha = alphalist[j]
  tau = knockoff.threshold(W3,fdr = alpha,offset = 1)
  rej3 = as.numeric(which(W3>=tau))
  power3[j] = power_cal(rej3, rand)
  fdr3[j] = fdr_cal(rej3, rand)
  
}



################## 4. GhostKnockoff ##################

### generate the summary statistics

ss_b = ss_se = ss_z = numeric(p)
for(l in 1:p){
  linearmodel <- lm(y ~ X[,l])
  ss_b[l] <- coef(summary(linearmodel))[2,1]
  ss_se[l] <- coef(summary(linearmodel))[2,2]
  ss_z[l] <- coef(summary(linearmodel))[2,3]
}

##### you can use the following code to replace the above one and generate summary stats using 
##### different sample sizes, for more details, see /different_sample_size/continuous.R
# linear_weights <- ((100:N.effect) - 99)^2
# N.all <- sample(100:N.effect, p, prob = linear_weights, replace = FALSE)
# N.median <- median(N.all)
# ss_b = ss_se = ss_z = numeric(p)
# for(l in 1:p){
#   linearmodel <- lm(y[1:N.all[l]] ~ X[1:N.all[l],l])
#   ss_b[l] <- coef(summary(linearmodel))[2,1]
#   ss_se[l] <- coef(summary(linearmodel))[2,2]
#   ss_z[l] <- coef(summary(linearmodel))[2,3]
# }



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
time4 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len){
  threshold = alphalist[j]
  rej4 <- which(GK.filter$q <= threshold)
  
  power4[j] = power_cal(rej4, rand)
  fdr4[j] = fdr_cal(rej4, rand)
}

############### 5. AnnoGK-simple ###############
start <- Sys.time()
set.seed(seed)
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

end <- Sys.time()
time5 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len){
  threshold = alphalist[j]
  rej5 <- which(GK.filter$q <= threshold)
  power5[j] = power_cal(rej5, rand)
  fdr5[j] = fdr_cal(rej5, rand)
}


############# 6. AnnoGK #############
start <- Sys.time()
set.seed(seed)
GK1ps_anno = GK_anno_M(Z, R, M, LD, N.effect)


beta <- GK1ps_anno$beta_final
T_0<-GK1ps_anno$T_0
T_k<-GK1ps_anno$T_k

GK.filter<-GhostKnockoff.filter(T_0,T_k)

end <- Sys.time()
time6 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len){
  threshold = alphalist[j]
  rej6 <- which(GK.filter$q <= threshold)
  power6[j] = power_cal(rej6, rand)
  fdr6[j] = fdr_cal(rej6, rand)
}



############# 7. AnnoGK with different sample sizes #############

start <- Sys.time()
set.seed(seed)
GK1ps_anno = GK_anno_dss(Z, R, M, LD, N.effect)

T_0_dss<-GK1ps_anno$T_0
T_k_dss<-GK1ps_anno$T_k

GK.filter<-GhostKnockoff.filter(T_0_dss,T_k_dss)

end <- Sys.time()

time7 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len){
  threshold = alphalist[j]
  rej7 <- which(GK.filter$q <= threshold)
  power7[j] = power_cal(rej7, rand)
  fdr7[j] = fdr_cal(rej7, rand)
}


############## 8. Ghosknockoff M = 5 ##############
start <- Sys.time()
M = 5 # number of knockoff copy
set.seed(seed)
fit.prelim <- GhostKnockoff.prelim(
  cor.G   = LD,
  M       = M,
  method  = "sdp" 
)
GK5_lasso <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method='lasso')
GK.filter<-GhostKnockoff.filter(GK5_lasso$T_0[[1]],GK5_lasso$T_k[[1]])

end <- Sys.time()
time8 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len){
  threshold = alphalist[j]
  rej8 <- which(GK.filter$q <= threshold)
  
  power8[j] = power_cal(rej8, rand)
  fdr8[j] = fdr_cal(rej8, rand)
}
 
############### 9. AnnoGK-simple M = 5 ###############
start <- Sys.time()
set.seed(seed)
GK5_M5_anno_ps  = GK_simple(Z = Z, 
                            R = R, 
                            M = M, 
                            LD = LD,
                            n = N.effect,
                            ts = 'lasso')
beta <- GK5_M5_anno_ps$beta
T_0<-abs(beta[1:p])
T_k<-abs(matrix(beta[-(1:p)],p,M))

GK.filter<-GhostKnockoff.filter(T_0,T_k)
end <- Sys.time()
time9 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len){
  threshold = alphalist[j]
  rej9 <- which(GK.filter$q <= threshold)
  power9[j] = power_cal(rej9, rand)
  fdr9[j] = fdr_cal(rej9, rand)
}


############# 10. AnnoGK M = 5 #############
start <- Sys.time()
set.seed(seed)
GK5ps_anno = GK_anno_M(Z, R, M, LD, N.effect)


beta <- GK5ps_anno$beta_final
T_0<-GK5ps_anno$T_0
T_k<-GK5ps_anno$T_k

GK.filter<-GhostKnockoff.filter(T_0,T_k)

end <- Sys.time()
time10 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len){
  threshold = alphalist[j]
  rej10 <- which(GK.filter$q <= threshold)
  power10[j] = power_cal(rej10, rand)
  fdr10[j] = fdr_cal(rej10, rand)
}


############# 11. AnnoGK different sample size M = 5 #############
start <- Sys.time()
set.seed(seed)
GK5ps_anno = GK_anno_dss(Z, R, M, LD, N.effect)


beta <- GK5ps_anno$beta_final
T_0<-GK5ps_anno$T_0
T_k<-GK5ps_anno$T_k

GK.filter<-GhostKnockoff.filter(T_0,T_k)

end <- Sys.time()
time11 <- as.numeric(difftime(end, start, units = "secs"))

for (j in 1:len){
  threshold = alphalist[j]
  rej11 <- which(GK.filter$q <= threshold)
  power11[j] = power_cal(rej11, rand)
  fdr11[j] = fdr_cal(rej11, rand)
}

