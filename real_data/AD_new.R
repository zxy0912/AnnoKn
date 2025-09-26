args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)

chrid = as.numeric(args[1])

source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")
packageVersion("Matrix") 


library(tidyr)
library(dplyr)
library(genio)
library(bigstatsr)
library(dendextend)

bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


################ load the summary statistics
ancestry0 = 'EUR'
ancestry1 = c('AFR')
M = 1
seed = 12345


ancestry = ancestry0
if(ancestry0 == 'AFR'){
  ancestry = 'AA'
}
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD/", ancestry,".txt")
sum_s <- read.table(path, header = TRUE, sep = '\t')
colnames(sum_s) <- c("SNPID", "SNP", "Chromsome","Position", "EffectAllele", "NonEffectAllele", "EffectAF", "Beta", "SE","Pval","Neff")
sum_s[[ancestry]] <- sum_s$Pval
print(length(unique(sum_s$SNP)) == nrow(sum_s))



