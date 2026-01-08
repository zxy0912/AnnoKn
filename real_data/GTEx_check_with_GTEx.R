library(arrow)
library(data.table)

allcelltype = c('Muscle_Skeletal', 'Whole_Blood', "Skin_Sun_Exposed_Lower_leg", "Artery_Tibial", 
                'Adipose_Subcutaneous', "Thyroid", "Nerve_Tibial", 
                "Skin_Not_Sun_Exposed_Suprapubic", "Lung", "Esophagus_Mucosa")

celltype = "Muscle_Skeletal"

sum_total <- numeric()
fdr_total <- numeric()
gene_total <- numeric()

for(celltype in allcelltype){
  print(celltype)
  
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/GTEx/GTEx_Analysis_v10_eQTL_updated/", 
                celltype, ".v10.eQTLs.signif_pairs.parquet")
  
  eqtls <- read_parquet(path)
  
  
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/GTEx/GTEx_Analysis_v10_eQTL_updated/", 
                celltype, ".v10.eGenes.txt.gz")
  
  data <- fread(path)
  
  data <- data[data$biotype == "protein_coding",]
  print(table(data$qval < 0.05))
  gene_total <- append(gene_total, nrow(data))
  
  sum_total <- append(sum_total, sum(data$qval < 0.05))
  
  
  eqtls$fdr_nominal <- p.adjust(eqtls$pval_nominal, method = "BH")
  print(table(eqtls$fdr_nominal < 0.1))
  eqtls <- eqtls[eqtls$fdr_nominal < 0.1,]
  eqtls <- eqtls[eqtls$gene_id %in% data$gene_id,]
  length(unique(eqtls$gene_id))
  fdr_total <- append(fdr_total, length(unique(eqtls$gene_id)))
  
}






allcelltype <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Muscle_Skeletal", "Skin_Sun_Exposed_Lower_leg", 
                 "Skin_Not_Sun_Exposed_Suprapubic", "Esophagus_Mucosa", "Colon_Sigmoid", "Colon_Transverse", "Small_Intestine_Terminal_Ileum")

for(celltype in allcelltype){
  print(celltype)
  
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/GTEx/GTEx_Analysis_v10_eQTL_updated/", 
                celltype, ".v10.eQTLs.signif_pairs.parquet")
  
  eqtls <- read_parquet(path)
  eqtls_df <- as.data.frame(eqtls)
  
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/GTEx/GTEx_Analysis_v10_eQTL_updated/", 
                celltype, ".v10.eQTLs.signif_pairs.csv")
  write.csv(eqtls_df, file = path)
  
}


