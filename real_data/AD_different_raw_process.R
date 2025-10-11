
library(dplyr)
GWAS1 <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_different/NG00075_Kunkle_etal_Stage1_P-val_only_results.txt",
                    header = TRUE)

head(GWAS1)
GWAS1$Ncase <- 21982
GWAS1$Ncontrol <- 41944

valid_nucleotides <- c("A", "C", "G", "T")
table(GWAS1$Effect_allele %in% valid_nucleotides & GWAS1$Non_Effect_allele %in% valid_nucleotides)
GWAS1_cleaned <- GWAS1 %>%
  filter(
    Effect_allele %in% valid_nucleotides &
      Non_Effect_allele %in% valid_nucleotides
  )
head(GWAS1_cleaned)
GWAS1_cleaned$N = 4/(1/GWAS1_cleaned$Ncase + 1/GWAS1_cleaned$Ncontrol)
# check is the hg19 base

write.table(GWAS1_cleaned,
            file = "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_different/IGAP.txt",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)



GWAS2 <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_different/3_UKB_AD_parental_meta_summary_output_June2019.txt",
                    header = TRUE)

head(GWAS2)
GWAS2$Ncase <- 27696 + 14338
GWAS2$Ncontrol <- 314278 - (27696 + 14338)

valid_nucleotides <- c("A", "C", "G", "T")
table(GWAS2$A1 %in% valid_nucleotides & GWAS2$A2 %in% valid_nucleotides)
GWAS2_cleaned <- GWAS2 %>%
  filter(
    A1 %in% valid_nucleotides &
      A2 %in% valid_nucleotides
  )
head(GWAS2_cleaned)
GWAS2_cleaned$N = 4/(1/GWAS2_cleaned$Ncase + 1/GWAS2_cleaned$Ncontrol)
bim <- read.table("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/simu_sep/1000G/EUR/1000G_EUR_maf1.bim")
# 
# head(GWAS2_cleaned)
# 
# common_snps <- merge(
#   bim, GWAS2_cleaned,
#   by.x = c("V1", "V4"),
#   by.y = c("CHR", "BP")
# )
# 
# n_common <- nrow(common_snps)
# print(n_common)
# check is the hg19 base

# bim_sub <- bim[, c("V1", "V4", "V2")]
# colnames(bim_sub) <- c("CHR", "BP", "rsid_bim")
# 
# GWAS2_with_rsid <- merge(
#   GWAS2_cleaned,
#   bim_sub,
#   by = c("CHR", "BP"),
#   all.x = TRUE
# )

sum(grepl("^rs", GWAS2_cleaned$SNP))

GWAS2_cleaned <- GWAS2_cleaned[, c(-9)]


# dup_snps
# [1] "19_404337"  "6_29838099" "6_32554197"
GWAS2_cleaned[498,]
GWAS2_cleaned <- GWAS2_cleaned[c(-498, -499, -963, -1008, -1009),]

write.table(GWAS2_cleaned,
            file = "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/AD_different/UKBB.txt",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

