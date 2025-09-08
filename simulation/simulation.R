

source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")


i = 11
seed <- 100*i

N.effect = 1000
p = 300 ## 200 500 1000
n0 = 30
amplitude = 3 ## 1 2 4
rho = 0.5
alphalist  <- seq(0.3, 0.05, by = -0.05)
fdr_seq  <- seq(0.3, 0.05, by = -0.05)
len = length(fdr_seq)


set.seed(seed)
power1 = numeric(len) ## knockoff
power2 = numeric(len) ## knockoff simple
power3 = numeric(len) ## GK M=1 pseudo sum
power4 = numeric(len) ## Ghostknockoff M=1 lasso min
power5 = numeric(len) ## Ghostknockoff M=5 pseudo sum
power6 = numeric(len) ## Ghostknockoff M=5 lasso min
power7 = numeric(len) ## GK M=1 pseudo sum with annotation simple
power8 = numeric(len) ## GK M=1 lasso min with annotation simple
power9 = numeric(len) ## GK M=1 pseudo sum with annotation
power10 = numeric(len) ## GK M=1 lasso min with annotation

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
  power1[j] = power_cal(rej1, rand)
  fdr1[j] = fdr_cal(rej1, rand)
  
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
  power2[j] = power_cal(rej2, rand)
  fdr2[j] = fdr_cal(rej2, rand)
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
  power3[j] = power_cal(rej3, rand)
  fdr3[j] = fdr_cal(rej3, rand)
  
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
  
  power4[j] = power_cal(rej4, rand)
  fdr4[j] = fdr_cal(rej4, rand)
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
  power7[j] = power_cal(rej7, rand)
  fdr7[j] = fdr_cal(rej7, rand)
}





############# 9. GK1 ps with annotation

GK1ps_anno = GK_anno(Z, R, M, LD, N.effect)


beta <- GK1ps_anno$beta_final
T_0<-abs(beta[1:p])
T_k<-abs(matrix(beta[-(1:p)],p,M))

GK.filter<-GhostKnockoff.filter(T_0,T_k)

for (j in 1:len){
  threshold = alphalist[j]
  rej9 <- which(GK.filter$q <= threshold)
  power9[j] = power_cal(rej9, rand)
  fdr9[j] = fdr_cal(rej9, rand)
}

