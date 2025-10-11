library(dplyr)

EUR_data <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/BC_published/selected_EUR_gwas.txt",
                       header = TRUE)

EUR_data$Ncase = 46785
EUR_data$Ncontrol = 42892

valid_nucleotides <- c("A", "C", "G", "T")
table(EUR_data$a0 %in% valid_nucleotides & EUR_data$a1 %in% valid_nucleotides)
EUR_data_cleaned <- EUR_data %>%
  filter(
    a0 %in% valid_nucleotides &
      a1 %in% valid_nucleotides
  )
head(EUR_data_cleaned)

EUR_data_cleaned$N = 4/(1/EUR_data_cleaned$Ncase + 1/EUR_data_cleaned$Ncontrol)


EUR_data_cleaned$phase3_1kg_id <- ifelse(
  EUR_data_cleaned$phase3_1kg_id == "NULL",
  NA,
  sub(":.*", "", EUR_data_cleaned$phase3_1kg_id)
)

table(is.na(EUR_data_cleaned$phase3_1kg_id))
EUR_data_cleaned <- EUR_data_cleaned[!is.na(EUR_data_cleaned$phase3_1kg_id),]

write.table(EUR_data_cleaned,
            file = "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/BC_published/EUR.txt",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)




EAS_data <- read.table("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/BC_published/selected_EAS_gwas.txt",
                       header = TRUE)

EAS_data$Ncase = 6269
EAS_data$Ncontrol = 6624

valid_nucleotides <- c("A", "C", "G", "T")
table(EAS_data$a0 %in% valid_nucleotides & EAS_data$a1 %in% valid_nucleotides)
EAS_data_cleaned <- EAS_data %>%
  filter(
    a0 %in% valid_nucleotides &
      a1 %in% valid_nucleotides
  )
head(EAS_data_cleaned)
EAS_data_cleaned$N = 4/(1/EAS_data_cleaned$Ncase + 1/EAS_data_cleaned$Ncontrol)

EAS_data_cleaned$phase3_1kg_id <- ifelse(
  EAS_data_cleaned$phase3_1kg_id == "NULL",
  NA,
  sub(":.*", "", EAS_data_cleaned$phase3_1kg_id)
)
table(is.na(EAS_data_cleaned$phase3_1kg_id))
EAS_data_cleaned <- EAS_data_cleaned[!is.na(EAS_data_cleaned$phase3_1kg_id),]

write.table(EAS_data_cleaned,
            file = "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GK_anno/BC_published/EAS.txt",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

