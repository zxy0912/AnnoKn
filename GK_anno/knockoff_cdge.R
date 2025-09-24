

library(KnockoffCDGE)
library(JuliaCall)
library(knockoffsr)
julia <- julia_setup("/home/xz527/julia-1.11.6/bin",installJulia = FALSE)
# ko <-knockoff_setup()

data("DGE_result")
data("Sigma")

head(DGE_result)

Z<-Z_calculation(DGE_result$pvalue,DGE_result$log2FoldChange)

set.seed(433)
CDGE_result<-KnockoffCDGE(Z,Sigma,M=5,n=23010,method="ME",Gene_info=DGE_result,verbose=TRUE,tol=0.001)
CDGE_result[1:20,]