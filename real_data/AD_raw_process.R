library(data.table)
library(dplyr)

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


AD_new_clean <- AD_new %>%
  group_by(SNP) %>%
  slice_min(order_by = p, n = 1, with_ties = FALSE) %>%
  ungroup()

AD_new_sorted <- arrange(AD_new_clean, chr, pos)



######################################### save the data ########################################
write.table(AD_new_sorted,"/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/GWAS2.txt", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)

path <- c("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/GWAS2.txt")
AD_new_sorted <- read.table(path, header = TRUE)

valid_nucleotides <- c("A", "C", "G", "T")
AD_new_cleaned <- AD_new_sorted %>%
  filter(
    A1 %in% valid_nucleotides &
      A2 %in% valid_nucleotides
  )
head(AD_new_cleaned)
### all SNPs have good, thus don't need to do this step
################################################################################################



################################################ GWAS 1 ####################################
AD_hg38 <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/GCST90027158_buildGRCh38.tsv",
                      sep = '\t', header = TRUE)

AD_hg38$SNP = paste(AD_hg38$chromosome, AD_hg38$base_pair_location, sep = "_")
AD_hg38$z = AD_hg38$beta/AD_hg38$standard_error
AD_hg38_recalculated <- AD_hg38 %>%
  mutate(
    # Calculate Z-score from beta and standard error
    z_score_calculated = beta / standard_error,
    
    # Calculate two-sided P-value from the new Z-score
    p_value_calculated = 2 * pnorm(abs(z_score_calculated), lower.tail = FALSE)
  )
AD_hg38_recalculated$N = AD_hg38_recalculated$n_cases + AD_hg38_recalculated$n_controls
AD_hg38_recalculated <- AD_hg38_recalculated[, c(18, 3, 4, 5, 6, 7, 11, 12, 20, 21, 22, 13, 1)]
colnames(AD_hg38_recalculated) <- c("SNP", "chr", "pos","A1", "A2", "maf", "beta", "se", "z", 
                                    "p", "N", "cases", "rs_id")


valid_nucleotides <- c("A", "C", "G", "T")
AD_hg38_clean <- AD_hg38_recalculated %>%
  filter(
    A1 %in% valid_nucleotides &
      A2 %in% valid_nucleotides
  )
head(AD_hg38_clean)

length(unique(AD_hg38_clean$SNP))

AD_hg38_final <- AD_hg38_clean %>%
  group_by(SNP) %>%
  slice_min(order_by = p, n = 1, with_ties = FALSE) %>%
  ungroup()

write.table(AD_hg38_final,"/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/GWAS1.txt", 
            row.names = FALSE, quote = FALSE, col.names = TRUE)


length(intersect(AD_hg38_recalculated$SNP, AD_new_sorted$SNP))




################################# divide into different chromosomes: #################################

study1 = 'GWAS1'
study2 = 'GWAS2'
M = 1

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/", study1,".txt")
sum_s <- read.table(path, header = TRUE)
sum_s <- arrange(sum_s, chr, pos)
print(length(unique(sum_s$SNP)) == nrow(sum_s))

for(chrid in 1:22){
  print(chrid)
  sum_s_chr = sum_s[sum_s$chr == chrid,]
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/bychr/", 
                study1, "_chr", chrid, ".txt")
  write.table(sum_s_chr, path, 
              row.names = FALSE, quote = FALSE, col.names = TRUE)
}



path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/", study2,".txt")
sum_s <- read.table(path, header = TRUE)
sum_s <- arrange(sum_s, chr, pos)
print(length(unique(sum_s$SNP)) == nrow(sum_s))

for(chrid in 1:22){
  print(chrid)
  sum_s_chr = sum_s[sum_s$chr == chrid,]
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/bychr/", 
                study2, "_chr", chrid, ".txt")
  write.table(sum_s_chr, path, 
              row.names = FALSE, quote = FALSE, col.names = TRUE)
}


##################### raw files for ldsc ##############################

study1 = 'GWAS1'
study2 = 'GWAS2'

library(dplyr)

all_file <- numeric()
for(chrid in 1:22){
  print(chrid)
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/bychr/", 
                study1, "_chr", chrid, ".txt")
  sum_s_chr <- read.table(path, header = TRUE)
  
  path <- paste("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/1KG_panel/hg38/EUR/bychr/1KG_chr", chrid,sep='')
  ref_chr <- read_plink(path)
  
  temp = sum_s_chr[, c("pos", "z", "A1","A2","N", "p", "chr")]
  temp_anno <- temp %>%
    left_join(ref_chr$bim[, c("pos", "id")], by = "pos")
  temp_anno <- temp_anno[, c(8,2,3,4,5,6,7)]
  temp_anno <- temp_anno[!is.na(temp_anno$id), ]
  all_file <- rbind(all_file, temp_anno)
  
}

colnames(all_file) <- c("SNP", "Z","A1","A2","N","Pval", "Chr")

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/LDSC/data/GWAS1_ldsc.txt")
write.table(all_file, path, 
            row.names = FALSE, quote = FALSE, col.names = TRUE)






all_file <- numeric()
for(chrid in 1:22){
  print(chrid)
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_new/bychr/", 
                study2, "_chr", chrid, ".txt")
  sum_s_chr <- read.table(path, header = TRUE)
  
  path <- paste("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/1KG_panel/hg38/EUR/bychr/1KG_chr", chrid,sep='')
  ref_chr <- read_plink(path)
  
  temp = sum_s_chr[, c("pos", "z", "A1","A2","N", "p", "chr")]
  temp_anno <- temp %>%
    left_join(ref_chr$bim[, c("pos", "id")], by = "pos")
  temp_anno <- temp_anno[, c(8,2,3,4,5,6,7)]
  temp_anno <- temp_anno[!is.na(temp_anno$id), ]
  all_file <- rbind(all_file, temp_anno)
  
}

colnames(all_file) <- c("SNP", "Z","A1","A2","N","Pval", "Chr")

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/LDSC/data/GWAS2_ldsc.txt")
write.table(all_file, path, 
            row.names = FALSE, quote = FALSE, col.names = TRUE)


