library(data.table)

AD_hg37 <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/PGCALZ2sumstatsExcluding23andMe.txt",
                      header = TRUE)

AD_hg37$SNP <- paste(AD_hg37$chr,":", AD_hg37$PosGRCh37, sep="")

liftover_bed = data.table(CHR = paste0("chr",AD_hg37$chr), Start = AD_hg37$PosGRCh37 - 1, 
                          End = AD_hg37$PosGRCh37, id = AD_hg37$SNP)
liftover_bed$Start <- as.integer(liftover_bed$Start)
liftover_bed$End <- as.integer(liftover_bed$End)
write.table(liftover_bed,"/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/lift.bed", 
            row.names = F, quote = F, col.names = F)


###### /gpfs/gibbs/pi/zhao/lx94/Software/liftOver lift.bed 
###### /gpfs/gibbs/pi/zhao/xz527/TWAS_fm/real_data/t2d/hg19ToHg38.over.chain 
###### lift_hg38.txt lift_exclude.txt


hg19to38 <- fread("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/lift_hg38.txt")
dim(hg19to38)

commonchr <- paste("chr",1:22,sep="")

colnames(hg19to38) <- c("CHR_38","start_38","end_38","SNP")

table(hg19to38$CHR_38 %in% commonchr)
index <- which(!hg19to38$CHR_38 %in% commonchr)
hg19to38 <- hg19to38[-index,]

AD_new <- merge(AD_hg37, hg19to38,by="SNP", sort=FALSE, all=TRUE)

summary(AD_new$end_38/AD_new$PosGRCh37)


index <- which(!is.na(AD_new$end_38))

AD_new <- AD_new[index,]

######################## remove bad SNPs and repeated SNPs by keeping the one with the smallest p-value
index = which(paste0("chr", AD_new$chr) != AD_new$CHR_38)
AD_new = AD_new[-index,]

AD_new$pos_id = paste(AD_new$chr, "_",AD_new$end_38,sep="")
dim(AD_new)
head(AD_new)
AD_new = AD_new[, c(12,2,10,11,4,5,6,7,8,3)]
colnames(AD_new) <- c("SNP", "chr", "start", "pos","A1", "A2", "z", "p", "N", "pos_hg37")

length(unique(AD_new$SNP))
dim(AD_new)
library(dplyr)

AD_new_clean <- AD_new %>%
  group_by(SNP) %>%
  slice_min(order_by = p, n = 1, with_ties = FALSE) %>%
  ungroup()

######################################### save the data ########################################
write.table(AD_new_clean,"/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/GWAS2.txt", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)
################################################################################################



################################################ GWAS 1 ####################################
AD_hg38 <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/GCST90027158_buildGRCh38.tsv",
                      sep = '\t', header = TRUE)

AD_hg38$SNP = paste(AD_hg38$chromosome, AD_hg38$base_pair_location)
AD_hg38 <- AD_hg38[, c(18, 3, 4, 5, 6, 7, 11, 12, 13, 14)]
