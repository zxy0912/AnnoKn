

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
ancestry = 'LAT'

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/SCZ/", ancestry, "_summ.tsv") 
sum_s <- read.table(path, sep="\t", header=TRUE, stringsAsFactors = FALSE)

sum_s <- sum_s[sum_s$PVAL<0.001,]

# path = paste("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/real_data/ancestry/",ancestry,"/T2D2024/t2d_ss_2024_sub", sep="")
# write.table(t2d_ss_2024, path, row.names = F, col.names = T, quote = F)

library("qqman")

p <- manhattan(sum_s, chr = "CHROM", bp = "POS", p = "PVAL", snp = "ID", main=paste("Manhattan plot for", ancestry, "ancestry"), 
               col = c("red4","red","orange","gold","darkgreen","forestgreen","yellowgreen","darkcyan","darkblue","royalblue1", "blue","dodgerblue","deepskyblue","skyblue1","purple4","darkmagenta","violetred","hotpink","palevioletred","lightpink","chocolate4","lightgray","gray28"), 
               ylim=c(0,42),chrlabs = NULL,suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08),highlight = NULL, 
               logp = TRUE, annotatePval = NULL, annotateTop = TRUE)

p





