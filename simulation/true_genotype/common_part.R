

bim1 <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/EUR_0.05maf/EUR_2W_new.bim")

bim2 <- read.table("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/simu_sep/1000G/EUR/1000G_EUR_maf1.bim")

common_snp <- intersect(bim1$V2, bim2$V2)

write.table(common_snp, 
            file = "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/EUR_0.05maf/EUR_common.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)





n = 10000

library(tidyr)
library(dplyr)
library(genio)

bedNA <- function(bed1){
  for(j in 1:ncol(bed1)){
    temp <- bed1[,j]
    temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
    bed1[,j] <- temp
    #print(j)
  }
  return(bed1)
}


maf_cal <- function(x){
  freq <- sum(x, na.rm = TRUE)/(2 * sum(!is.na(x)))
  maf <- ifelse(freq < 0.5, freq, 1 - freq)
  return(maf)
}

ancestry = 'EUR'
chrid = '1'
path <- paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/",ancestry,"_0.05maf/bychr/chr_",chrid,"_2w_common")
genotype <- read_plink(path)


# select the X from a LD block:
ld_block <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/EUR/fourier_ls-all.bed", header = TRUE)
ld_block = ld_block[ld_block$chr == 'chr1',]
i = 1

#################### find proper regions

snps <- numeric()
snp_clusters <- numeric()
for(i in 1: nrow(ld_block)){
  
  print(i)
  
  start = ld_block[i,2]
  end = ld_block[i,3]
  index = which(genotype$bim$pos >= start & genotype$bim$pos <= end)
  
  X = t(genotype$X[index, 1:n])
  maf <- apply(X, 2, maf_cal)
  X = bedNA(X)
  
  snps <- append(snps, ncol(X))
  
  #### use hierarchical clustering to cluster SNPs
  
  R = cor(X[1:1000,])
  R2 <- R^2
  d <- as.dist(1 - R2)  # higher distance = lower correlation
  hc <- hclust(d, method = "average")
  plot(hc, labels = FALSE, main = "SNP clustering by LD")
  clusters <- cutree(hc, h = 0.75)  # e.g., cluster SNPs with r² > 0.25
  length(table(clusters))
  
  snp_clusters = append(snp_clusters, length(table(clusters)))
}

snp_number <- data.frame(index = 1:nrow(ld_block),
                         snps = snps,
                         snp_clusters = snp_clusters)
write.table(snp_number, 
            "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/EUR_0.05maf/real_simulation/region_information/chr1_common.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)

good_iter <- which(snp_number$snp_clusters > 100)




######## find 2MB windows:



ancestry = 'EUR'
chrid = '1'
path <- paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/",ancestry,"_0.05maf/bychr/chr_",chrid,"_2w_common")
genotype <- read_plink(path)

start_chr = min(genotype$bim$pos)
end_chr = max(genotype$bim$pos)

width = 2000000
n = 10000

snps <- numeric()
snp_clusters <- numeric()
start = start_chr
start_region <- numeric()
end_region <- numeric()

while (start < end_chr) {
  
  print(start)
  end = min(start + width, end_chr)
  
  index = which(genotype$bim$pos >= start & genotype$bim$pos <= end)
  
  if(length(index) == 0){
    start = start + width
    next
  }
  
  start_region <- append(start_region, start)
  end_region <- append(end_region, end)
  
  X = t(genotype$X[index, 1:n])
  maf <- apply(X, 2, maf_cal)
  X = bedNA(X)
  
  snps <- append(snps, ncol(X))
  
  #### use hierarchical clustering to cluster SNPs
  
  R = cor(X[1:1000,])
  R2 <- R^2
  d <- as.dist(1 - R2)  # higher distance = lower correlation
  hc <- hclust(d, method = "average")
  plot(hc, labels = FALSE, main = "SNP clustering by LD")
  clusters <- cutree(hc, h = 0.75)  # e.g., cluster SNPs with r² > 0.25
  length(table(clusters))
  
  snp_clusters = append(snp_clusters, length(table(clusters)))
  
  start = start + width
}

snp_number <- data.frame(index = 1:length(snps),
                         start = start_region,
                         end = end_region,
                         snps = snps,
                         snp_clusters = snp_clusters)
write.table(snp_number, 
            "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/UKBB/EUR_0.05maf/real_simulation/region_information/chr1_common_2MB.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)


