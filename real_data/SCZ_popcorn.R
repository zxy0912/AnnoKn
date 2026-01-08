

library(data.table)
bim=fread("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/simu_sep/1000G/EAS/1000G_EAS_step3.bim")
bim1=bim[duplicated(bim$V2),]
bim2=bim[!bim$V2 %in% bim1$V2,]
c=bim2
c1=c[!((c$V5=="G")&(c$V6 =="C")|(c$V5=="C")&(c$V6 =="G")|(c$V5=="A")&(c$V6 =="T")|(c$V5=="T")&(c$V6 =="A")),]

write.table(c1$V2,"/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/simu_sep/1000G/EAS/filter_SNP.txt", row.names=F,col.names=F,quote=F)



library(dplyr)

# ancestry_all <- c("EUR", "AFR")

common_snp = c("A","T","G","C")

ancestry = 'EUR'
ancestry0 = 'european'
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/PGC3_SCZ_wave3.", ancestry0, ".autosome.public.v3.vcf.tsv") 
sum_s <- read.table(path,
                    sep="\t", header=TRUE, comment.char="#", stringsAsFactors=FALSE)


path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/EUR_summ.tsv") 
write.table(sum_s, file = path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

colnames(sum_s) <- c("chr", "rsid", "pos", "a1", "a2", "af_case", "af_control", "info",
                     "beta", "SE","Pval","N_case", "N_control", "N")
table(sum_s$a1 %in% common_snp & sum_s$a2 %in% common_snp)
sum_s$AF <- 1 - (sum_s$af_case*sum_s$N_case + sum_s$af_control*sum_s$N_control)/(sum_s$N_case + sum_s$N_control)
# sum_s <- sum_s[, c(2, 5, 6, 7, 11, 8, 9)]
output_path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/Popcorn/data_SCZ/", ancestry,"summs.txt")
sum_s_formatted <- sum_s %>%
  dplyr::select(
    rsid = rsid,
    a1 = a1,
    a2 = a2,
    AF = AF,
    N = N,
    beta = beta,
    SE = SE
  )
readr::write_tsv(sum_s_formatted, output_path)




ancestry = 'EAS'
ancestry0 = 'asian'
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/PGC3_SCZ_wave3.", ancestry0, ".autosome.public.v3.vcf.tsv") 
sum_s <- read.table(path,
                    sep="\t", header=TRUE, comment.char="#", stringsAsFactors=FALSE)

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/EAS_summ.tsv") 
write.table(sum_s, file = path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

colnames(sum_s) <- c("chr", "rsid", "pos", "a1", "a2", "af_case", "af_control", "info",
                     "beta", "SE","Pval","N_case", "N_control", "N")
table(sum_s$a1 %in% common_snp & sum_s$a2 %in% common_snp)
sum_s$AF <- 1 - (sum_s$af_case*sum_s$N_case + sum_s$af_control*sum_s$N_control)/(sum_s$N_case + sum_s$N_control)
# sum_s <- sum_s[, c(2, 5, 6, 7, 11, 8, 9)]
output_path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/Popcorn/data_SCZ/", ancestry,"summs.txt")
sum_s_formatted <- sum_s %>%
  dplyr::select(
    rsid = rsid,
    a1 = a1,
    a2 = a2,
    AF = AF,
    N = N,
    beta = beta,
    SE = SE
  )
readr::write_tsv(sum_s_formatted, output_path)




ancestry = 'AFR'
ancestry0 = 'afram'
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/PGC3_SCZ_wave3.", ancestry0, ".autosome.public.v3.vcf.tsv") 
sum_s <- read.table(path,
                    sep="\t", header=TRUE, comment.char="#", stringsAsFactors=FALSE)
sum_s$N_case <- 6152
sum_s$N_control <- 3918
sum_s$N <- (6152 + 3918)/2
# number obtained from https://www.nature.com/articles/s41380-019-0517-y

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/AFR_summ.tsv") 
write.table(sum_s, file = path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

colnames(sum_s) <- c("chr", "rsid", "pos", "a1", "a2", "af_case", "af_control", "info",
                     "beta", "SE","Pval","N_case", "N_control", "N")
table(sum_s$a1 %in% common_snp & sum_s$a2 %in% common_snp)
sum_s$AF <- 1 - (sum_s$af_case*sum_s$N_case + sum_s$af_control*sum_s$N_control)/(sum_s$N_case + sum_s$N_control)
# sum_s <- sum_s[, c(2, 5, 6, 7, 11, 8, 9)]
output_path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/Popcorn/data_SCZ/", ancestry,"summs.txt")
sum_s_formatted <- sum_s %>%
  dplyr::select(
    rsid = rsid,
    a1 = a1,
    a2 = a2,
    AF = AF,
    N = N,
    beta = beta,
    SE = SE
  )

readr::write_tsv(sum_s_formatted, output_path)




ancestry = 'LAT'
ancestry0 = 'latino'
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/PGC3_SCZ_wave3.", ancestry0, ".autosome.public.v3.vcf.tsv") 
sum_s <- read.table(path,
                    sep="\t", header=TRUE, comment.char="#", stringsAsFactors=FALSE)
sum_s$N_case <- 1234
sum_s$N_control <- 3090
sum_s$N <- (1234 + 3090)/2

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/LAT_summ.tsv") 
write.table(sum_s, file = path, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



########### make the manhatten plot
ancestry = 'EUR'

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/", ancestry, "_summ.tsv") 
sum_s <- read.table(path, sep="\t", header=TRUE, stringsAsFactors = FALSE)

sum_s <- sum_s[sum_s$PVAL<0.001,]

# path = paste("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/real_data/ancestry/",ancestry,"/T2D2024/t2d_ss_2024_sub", sep="")
# write.table(t2d_ss_2024, path, row.names = F, col.names = T, quote = F)

library("qqman")

p <- manhattan(sum_s, chr = "CHROM", bp = "POS", p = "PVAL", snp = "ID", main="", #paste("Manhattan plot for", ancestry, "ancestry"), 
               col = c("red4","red","orange","gold","darkgreen","forestgreen","yellowgreen","darkcyan","darkblue","royalblue1", "blue","dodgerblue","deepskyblue","skyblue1","purple4","darkmagenta","violetred","hotpink","palevioletred","lightpink","chocolate4","lightgray","gray28"), 
               ylim=c(0,42),chrlabs = NULL,suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08),highlight = NULL, 
               logp = TRUE, annotatePval = NULL, annotateTop = TRUE)

p










########## check the number of risk regions:




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
ancestry1 = c('EAS',"AFR","LAT")


ancestry = ancestry0
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/", ancestry, "_summ.tsv")
sum_s <- read.table(path, sep="\t", header=TRUE, stringsAsFactors = FALSE)
# we also checked that all SNPs are in (A,T,G,C)
colnames(sum_s) <- c("Chromsome","SNPID", "Position", "EffectAllele", "NonEffectAllele", "NonEffectAF_case", 
                     "NonEffectAF_control", "INFO", "Beta", "SE","Pval","N_case","N_control","Neff")
sum_s[[ancestry]] <- sum_s$Pval
print(length(unique(sum_s$SNPID)) == nrow(sum_s))

# ,'HIS','SAS','EAS','AFA'
for(ancestry in ancestry1){
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/", ancestry, "_summ.tsv")
  temp <- read.table(path, sep="\t", header=TRUE, stringsAsFactors = FALSE)
  colnames(temp) <- c("Chromsome","SNPID", "Position", "EffectAllele", "NonEffectAllele", "NonEffectAF_case", 
                      "NonEffectAF_control", "INFO", "Beta", "SE","Pval","N_case","N_control","Neff")
  
  message(ancestry, ": unique SNP? ", length(unique(temp$SNPID)) == nrow(temp))
  index <- match(sum_s$SNP, temp$SNP)
  print(table(is.na(index)))
  print(table(temp$SNP[index] == sum_s$SNP))
  message(ancestry, ": matched = ", sum(!is.na(index)), "/", length(index))
  sum_s[[ancestry]] <- temp$Pval[index]
}



number_all <- numeric()
risk_region_all1 <- c()

for(chrid in 1:22){
  
  print(paste("chr:", chrid))
  
  sums_chr = sum_s[sum_s$Chromsome == chrid,]
  sums_chr <- sums_chr[order(sums_chr$Position),]
  
  ################ load the reference panel from 1KG
  
  # ancestry = 'EUR'
  if(ancestry0 == 'EUR'){
    path <- paste("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/simu_sep/1000G/AFR/bychr_", ancestry0, "/1KG_chr", chrid,sep='')
  }else if(ancestry0 == 'EAS'){
    path <- paste("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/simu_sep/1000G/EAS/bychr_", ancestry0, "/1KG_chr", chrid,sep='')
  }
  ref_chr <- read_plink(path)
  
  ############### find the common part of three datasets:
  
  common_id <- intersect(ref_chr$bim$pos, sums_chr$Position)
  index1 = match(common_id, ref_chr$bim$pos)
  index2 = match(common_id, sums_chr$Position)
  
  bim = ref_chr$bim[index1, ]
  X = ref_chr$X[index1, ]
  table(rownames(X) == bim$id)
  sums_chr = sums_chr[index2, ]
  table(bim$pos == sums_chr$Position)
  
  
  ################### define risk regions
  
  len_half = 250000
  # risk_id <- which(sums_chr$Pval < 5 * 10^-8)
  risk_id <- which(sums_chr$EUR < 5 * 10^-8)
  pos_id = unlist(sums_chr$Position[risk_id])
  pos_diff <- diff(pos_id)
  where_sep = as.numeric(c(0, which(pos_diff > 2 * len_half)))
  risk_region = numeric()
  
  
  if(length(where_sep) > 1){
    for(i in 1:(length(where_sep)-1)){
      risk_region <- rbind(risk_region, c(pos_id[where_sep[i] + 1] - len_half, pos_id[where_sep[i+1]] + len_half))
    }
    risk_region <- rbind(risk_region, c(pos_id[where_sep[length(where_sep)] + 1] - len_half, pos_id[length(pos_id)] + len_half))
  }else if(length(where_sep) == 1){
    risk_region <- rbind(risk_region, c(pos_id[where_sep[length(where_sep)] + 1] - len_half, pos_id[length(pos_id)] + len_half))
  }
  
  print(risk_region)
  risk_region <- as.data.frame(risk_region)
  
  if(ncol(risk_region) > 1){
    risk_region$chr = chrid
    risk_region_all1 <- rbind(risk_region_all1, risk_region)
  }
}

