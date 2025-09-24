
library(dplyr)

ancestry_all <- c("EUR", "AFR")

for(ancestry in ancestry_all){
  if(ancestry == 'AFR'){
    ancestry0 = 'AA'
  }else{
    ancestry0 = ancestry
  }
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/Height/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_", ancestry0) 
  sum_s <- read.table(path, header = TRUE, sep = '\t')
  colnames(sum_s) <- c("SNPID", "rsid", "chr","pos", "a1", "a2", "af", "beta", "SE","Pval","N")
  sum_s <- sum_s[, c(2, 5, 6, 7, 11, 8, 9)]
  output_path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/Popcorn/data/", ancestry,"summs.txt")
  sum_s_formatted <- sum_s %>%
    dplyr::select(
      rsid = rsid,
      a1 = a1,
      a2 = a2,
      af = af,
      N = N,
      beta = beta,
      SE = SE
    )
  
  readr::write_tsv(sum_s_formatted, output_path)
}



########### make the manhatten plot
ancestry = 'AFR'
if(ancestry == 'AFR'){
  ancestry0 = 'AA'
}else{
  ancestry0 = ancestry
}
path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/Height/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_", ancestry0) 
sum_s <- read.table(path, header = TRUE, sep = '\t')

# t2d_ss_2024 <- t2d_ss_2024_all[t2d_ss_2024_all$Pval<0.0001,]

# path = paste("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/real_data/ancestry/",ancestry,"/T2D2024/t2d_ss_2024_sub", sep="")
# write.table(t2d_ss_2024, path, row.names = F, col.names = T, quote = F)

library("qqman")

p <- manhattan(sum_s, chr = "CHR", bp = "POS", p = "P", snp = "RSID", main=paste("Manhattan plot for",ancestry,"ancestry"), 
               col = c("red4","red","orange","gold","darkgreen","forestgreen","yellowgreen","darkcyan","darkblue","royalblue1", "blue","dodgerblue","deepskyblue","skyblue1","purple4","darkmagenta","violetred","hotpink","palevioletred","lightpink","chocolate4","lightgray","gray28"), 
               ylim=c(0,350),chrlabs = NULL,suggestiveline = -log10(5e-06), genomewideline = -log10(5e-08),highlight = NULL, 
               logp = TRUE, annotatePval = NULL, annotateTop = TRUE)

p





