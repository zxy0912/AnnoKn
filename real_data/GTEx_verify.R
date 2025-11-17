


knockoff.threshold <- function (W, fdr = 0.1, offset = 1){
  if (offset != 1 && offset != 0) {
    stop("Input offset must be either 0 or 1")
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t) (offset + sum(W <= -t))/max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}




# 

library(data.table)

sum_all_ct <- data.frame(matrix(nrow = 22, ncol = 0))
# number_gene_ct <- data.frame(matrix(nrow = 22, ncol = 0))
find1_all <- list()
find2_all <- list()


allcelltype = c('Muscle_Skeletal', 'Whole_Blood', "Skin_Sun_Exposed_Lower_leg", "Artery_Tibial", 
                'Adipose_Subcutaneous', "Thyroid", "Nerve_Tibial", 
                "Skin_Not_Sun_Exposed_Suprapubic", "Lung", "Esophagus_Mucosa")

# 

type = 1
alpha = 0.1
pc = 5
lnc = FALSE


annot_all = fread(paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/LDSC/data/GRCh38/baselineLD_v2.2/baselineLD.", 1, ".annot.gz"))

cats = c("SNP", "Promoter_UCSC", "Enhancer_Hoffman", "Coding_UCSC", "TFBS_ENCODE")
if(type == 1){
  cats = c("SNP", "Promoter_UCSC", "Enhancer_Hoffman", "TFBS_ENCODE")
}else if(type == 2){
  cats <- c(
    "SNP",
    # QTL
    "GTEx_eQTL_MaxCPP",
    "BLUEPRINT_H3K27acQTL_MaxCPP",
    "BLUEPRINT_H3K4me1QTL_MaxCPP",
    "BLUEPRINT_DNA_methylation_MaxCPP",
    
    # Promoter / TSS
    "Promoter_UCSC",
    "TSS_Hoffman",
    "PromoterFlanking_Hoffman",
    "Human_Promoter_Villar",
    "Ancient_Sequence_Age_Human_Promoter",
    "Human_Promoter_Villar_ExAC",
    
    # Enhancer / TFBS
    "Enhancer_Andersson",
    "Enhancer_Hoffman",
    "SuperEnhancer_Hnisz",
    "WeakEnhancer_Hoffman",
    "TFBS_ENCODE",
    "Human_Enhancer_Villar",
    "Ancient_Sequence_Age_Human_Enhancer",
    
    # Histone marks
    "H3K27ac_Hnisz",
    "H3K27ac_PGC2",
    "H3K4me1_peaks_Trynka",
    "H3K4me1_Trynka",
    "H3K4me3_peaks_Trynka",
    "H3K4me3_Trynka",
    "H3K9ac_peaks_Trynka",
    "H3K9ac_Trynka"
  )
}else if(type == 3){
  cats <- colnames(annot_all)[c(-1, -2, -4, -5)]
}

a_dim <- length(cats) - 1


lambdas_all <- list()
significant1_all <- c()
significant2_all <- c()
baselin_gene_all <- c()


for(celltype in allcelltype){
  
  significant1 <- c()
  significant2 <- c()
  baselin_gene <- c()
  
  lambda_s <- numeric()

  sum_all <- numeric()
  sum1_all = c()
  sum2_all = c()
  find1 <- c()
  find2 <- c()
  for(chrid in (1:22)){
    # print(chrid)
    sum1 = c()
    sum2 = c()
    path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/chr_", chrid, "_ct_", celltype, "_type_", type, ".RData")
    
    if (!file.exists(path)) {
      message("not exist: ", path)
      break
    }
    load(path)
    
    # print(length(all_result))
    
    sum_all <- append(sum_all, length(all_result))
    
    for(gene_i in 1:length(all_result)){
      # print(gene_i)
      adakn_result = all_result[[gene_i]]
      
      
      baselin_gene <-  rbind(baselin_gene, adakn_result$gene_info)
      
      W1 = adakn_result$ghostknockoff
      tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
      rej1 = as.numeric(which(W1>=tau))
      # print("rej1:")
      # print(rej1)
      if(length(rej1 > 0)){
        sum1 = append(sum1, gene_i)
        find1 <- append(find1, paste(chrid, gene_i))
        significant1 <- rbind(significant1, adakn_result$gene_info)
      }
      result_raw <- adakn_result$knockoff_anno
      
      
      if(length(unlist(result_raw$lambda_s)) < a_dim){
        print(paste("less than", a_dim,"annotation at celltype", celltype, "chr", adakn_result$chr, "gene", adakn_result$gene_i))
        # print(celltype)
        # print(adakn_result$chr)
        # print(adakn_result$gene_i)
        b = apply(adakn_result$R, 2, var)
        a = numeric(a_dim)
        a[b!=0] = unlist(result_raw$lambda_s)
        a[b==0] = 0
        lambda_s <- rbind(lambda_s, a)
      }else{
        lambda_s <- rbind(lambda_s, unlist(result_raw$lambda_s))
      }
      
      
      beta = result_raw$beta
      p = length(beta)/2
      T = abs(beta[1:p])
      T_tilde = abs(beta[(p+1):(2*p)])
      T_max = pmax(T,T_tilde)
      W = T - T_tilde
      # hist(W, breaks = 50)
      tau = knockoff.threshold(W,fdr = alpha,offset = 1)
      rej2 = as.numeric(which(W >= tau))
      
      # print("rej3:")
      # print(rej3)
      if(length(rej2 > 0)){
        sum2 = append(sum2, gene_i)
        find2 <- append(find2, paste(adakn_result$chr, adakn_result$gene_i))
        significant2 <- rbind(significant2, adakn_result$gene_info)
      }
    }
    sum1_all <- append(sum1_all, length(sum1))
    sum2_all <- append(sum2_all, length(sum2))
   
  }
  
  # cbind(number_gene, sum_all)
  
  print(length(find1))
  print(length(find2))
  
  sum_all_ct[[celltype]] = sum_all
  # number_gene_ct[[celltype]] = number_gene
  
  find1_all <- append(find1_all, list(find1))
  find2_all <- append(find2_all, list(find2))
  
  lambdas_all <- append(lambdas_all, list(lambda_s))

  significant1$celltype = celltype
  significant2$celltype = celltype
  baselin_gene$celltype = celltype
  
  significant1_all = rbind(significant1_all, significant1)
  significant2_all = rbind(significant2_all, significant2)
  baselin_gene_all = rbind(baselin_gene_all, baselin_gene)
  
}




significant1_all$hgnc_symbol[significant1_all$hgnc_symbol == ''] = significant1_all$ensembl_gene_id[significant1_all$hgnc_symbol == '']
significant2_all$hgnc_symbol[significant2_all$hgnc_symbol == ''] = significant2_all$ensembl_gene_id[significant2_all$hgnc_symbol == '']


improve_prop = (sapply(find2_all, length) - sapply(find1_all, length))/sapply(find1_all, length)

allcelltypename <- c("Muscle Skeletal", "Whole Blood", "Skin Sun Exposed (Lower leg)",
                     "Artery Tibial", "Adipose Subcutaneous", "Thyroid", "Nerve Tibial",
                     "Skin Not Sun Exposed (Suprapubic)", "Lung", "Esophagus Mucosa")

data.frame(tissue = allcelltypename, 
           knockoff = sapply(find1_all, length), 
           annokn = sapply(find2_all, length))


########### check the number of genes from different numbers of tissues ##############

library(dplyr)
gene_counts1 <- significant1_all %>%
  group_by(hgnc_symbol) %>%
  summarise(n_celltypes = n_distinct(celltype))
table(gene_counts1$n_celltypes)


gene_counts2 <- significant2_all %>%
  group_by(hgnc_symbol) %>%
  summarise(n_celltypes = n_distinct(celltype))
table(gene_counts2$n_celltypes)

sum(gene_counts1$n_celltypes > 1)
sum(gene_counts2$n_celltypes > 1)








########################### figure 1 plot #########################


library(ggplot2)
library(dplyr)

#data_plot <- data_plot %>%
#  arrange(desc(improve_prop))

allcelltypename <- c("Muscle Skeletal", "Whole Blood", "Skin Sun Exposed (Lower leg)",
                     "Artery Tibial", "Adipose Subcutaneous", "Thyroid", "Nerve Tibial",
                     "Skin Not Sun Exposed (Suprapubic)", "Lung", "Esophagus Mucosa")

improve_prop = (sapply(find2_all, length) - sapply(find1_all, length))/sapply(find1_all, length)
data_plot = data.frame(improve_prop = improve_prop,
                       allcelltype = allcelltypename)


data_plot$allcelltype <- factor(data_plot$allcelltype, levels = data_plot$allcelltype)

p <- ggplot(data_plot, aes(x = allcelltype, y = improve_prop)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    x = "Tissue",
    y = "Increased proportion",
    title = NULL
      # "Increased proportion of detected eGenes by AnnoKn over Knockoffs by tissue"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11) 
  ) + 
  geom_text(aes(label = round(improve_prop, 3)), vjust = -0.3, hjust = 0.4)

p


ggsave("/gpfs/gibbs/pi/zhao/xz527/annoKn_plots/gtex_type1_ahpla0.1.pdf", 
       p, 
       width = 9, height = 6, units = "in", 
       bg = "white", device = cairo_pdf)


data_plot_0.1 = data_plot # this is from using alpha = 0.1, run this code again
data_plot_0.2 = data_plot # this is from using alpha = 0.2, run this code again






data_plot_0.1$threshold <- "q = 0.1"
data_plot_0.2$threshold <- "q = 0.2"


combined_data <- rbind(data_plot_0.1, data_plot_0.2)

tissue_order <- data_plot_0.1$allcelltype
combined_data$allcelltype <- factor(combined_data$allcelltype, levels = tissue_order)


p <- ggplot(combined_data, aes(x = allcelltype, y = improve_prop, fill = threshold)) +
  
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  geom_text(
    aes(label = round(improve_prop, 3)), 
    position = position_dodge(width = 0.9), 
    vjust = -0.4, 
    size = 3      
  ) +
  
  labs(
    x = "Tissue",
    y = "Increased proportion",
    title = NULL,
    fill = "Threshold"
  ) +
  
  scale_fill_brewer(palette = "Paired") + 
  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "top" 
  )

p


ggsave("/gpfs/gibbs/pi/zhao/xz527/annoKn_plots/gtex_type3.pdf", 
       p, 
       width = 9, height = 6, units = "in", 
       bg = "white", device = cairo_pdf)












#### count the common number of genes from a list of tissues:


count_common_genes <- function(data, tissues, id_col = "ensembl_gene_id", tissue_col = "celltype", symbol_col = "hgnc_symbol") {
  
  missing_tissues <- setdiff(tissues, unique(data[[tissue_col]]))
  if (length(missing_tissues) > 0) {
    stop(paste("Tissue(s) not found in data:", paste(missing_tissues, collapse = ", ")))
  }
  
  common_genes <- unique(data[data[[tissue_col]] == tissues[1], id_col])
  
  for (t in tissues[-1]) {
    genes_t <- unique(data[data[[tissue_col]] == t, id_col])
    common_genes <- intersect(common_genes, genes_t)
  }
  
  n_common <- length(common_genes)
  cat("Number of common genes across", length(tissues), "tissues:", n_common, "\n")
  
  common_symbols <- unique(data[data[[id_col]] %in% common_genes, symbol_col])
  
  return(list(
    n_common = n_common,
    common_genes = common_genes,
    common_symbols = common_symbols
  ))
}


result_skin <- count_common_genes(
  significant2_all,
  tissues = c("Adipose_Visceral_Omentum", 'Adipose_Subcutaneous')
)
result_skin$n_common
head(result_skin$common_symbols)




################ double check the list of genes from a list of tissues.

genes_sun <- significant2_all$hgnc_symbol[significant2_all$celltype == "Skin_Sun_Exposed_Lower_leg"]
genes_not_sun <- significant2_all$hgnc_symbol[significant2_all$celltype == "Skin_Not_Sun_Exposed_Suprapubic"]
common_genes <- intersect(genes_sun, genes_not_sun)
length(common_genes)



################## check the enrichment with pathways: ##################
# library(devtools)
# install_github("wjawaid/enrichR")

library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs_tissue <- dbs[grep("tissue", dbs$libraryName, ignore.case = TRUE), ]
dbs <- c("Jensen_TISSUES", "ARCHS4_Tissues")
scan = 1


celltype = allcelltype[10]
print(celltype)
input = significant1_all$hgnc_symbol[significant1_all$celltype == celltype]
background = baselin_gene_all$hgnc_symbol[baselin_gene_all$celltype == celltype]
enriched2 <- enrichr(input, dbs, background = background)
print(head(enriched2[[dbs[scan]]][, c(1,3)], 20))
target_terms <- enriched2[[dbs[scan]]][grep("Esophagus", enriched2[[dbs[scan]]]$Term, ignore.case = TRUE), ]
target_terms[, c(1,3)]

input = significant2_all$hgnc_symbol[significant2_all$celltype == celltype]
enriched2 <- enrichr(input, dbs, background = background)
print(head(enriched2[[dbs[scan]]][, c(1,3)], 20))
target_terms <- enriched2[[dbs[scan]]][grep("Esophagus", enriched2[[dbs[scan]]]$Term, ignore.case = TRUE), ]
target_terms[, c(1,3)]


tissue_target <- c("Skeletal muscle", "Blood", "skin", "Coronary artery","Adipose tissue","Thyroid gland",
                   "Brain","skin","Lung","Esophagus")
allcelltype = c('Muscle_Skeletal', 'Whole_Blood',"Skin_Sun_Exposed_Lower_leg","Artery_Tibial", 
                'Adipose_Subcutaneous', "Thyroid", "Nerve_Tibial", 
                "Skin_Not_Sun_Exposed_Suprapubic", "Lung", "Esophagus_Mucosa")
allcelltypename <- c("Muscle Skeletal", "Whole Blood", "Skin Sun Exposed (Lower leg)",
                     "Artery Tibial", "Adipose Subcutaneous", "Thyroid", "Nerve Tibial",
                     "Skin Not Sun Exposed (Suprapubic)", "Lung", "Esophagus Mucosa")
result1 <- c(0.001102511, 0.0002752536, 0.1631501, 0.1414703, 0.0004044594, 5.399899e-06, 0.2345448, 0.08218199, 0.07763491, 0.01199245)
# type 1:
result2 <- c(0.0003043275, 0.0001758909, 0.02153168, 0.1687078, 7.180866e-05, 5.659479e-09, 0.005001548, 0.00424484, 0.006609615, 4.038481e-05)
# type 2:
result3 <- c(4.753923e-05, 4.229528e-05, 0.1208330, 0.06237147, 0.0002535993, 2.925358e-10, 0.03713819, 0.01093074, 0.01500972, 1.688129e-05)
# type 3:
result4 <- c(9.255167e-05, 1.904593e-05, 0.04613392, 0.06381061, 0.0004747116, 3.441058e-07, 0.01644105, 0.004317677, 0.001478324, 7.368574e-06)


logp1 <- -log10(result1)
logp3 <- -log10(result4)

library(ggplot2)
library(reshape2)

df <- data.frame(
  Tissue = factor(allcelltypename, levels = rev(allcelltypename)),
  Knockoff = logp1,
  AnnoKn = logp3
)

df_melt <- melt(df, id.vars = "Tissue", variable.name = "Type", value.name = "negLogP")


p <- ggplot(df_melt, aes(x = Tissue, y = negLogP, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(x = "", y = expression(-log[10](p-value)),
       title = NULL) +
  theme_minimal(base_size = 13) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  annotate("text", x = 7, y = -log10(0.05)+0.3, label = "p = 0.05", color = "red", hjust = 0)

ggsave("/gpfs/gibbs/pi/zhao/xz527/annoKn_plots/enrichr_type3_ahpla0.2.pdf", 
       p, 
       width = 10, height = 5, units = "in", 
       bg = "white", device = cairo_pdf)



################## check the lambdas: ##################


anno_name <- cats[-1]
anno_name <- c("Promoter", "Enhancer", "TFBS")
allcelltypename

lambda_result_small <- lambda_result_large <- data.frame(matrix(NA, 
                                         nrow = length(allcelltypename), 
                                         ncol = length(anno_name),
                                         dimnames = list(allcelltypename, anno_name)))

for(j in 1:3){
  
  for(i in 1:10){
    celltype = allcelltypename[i]
    name = anno_name[j]
    
    limited <- lambdas_all[[i]][, j]
    print(length(limited))
    limited <- limited[limited!=0]
    # limited <- limited[which(abs(limited) < 2.5)]
    print(length(limited))
    df <- data.frame(x = limited)
    
    print(sum(limited < -1))
    lambda_result_small[i,j] = sum(limited < -1)
    print(sum(limited > 1))
    lambda_result_large[i,j] = sum(limited > 1)
    
    p = ggplot(df, aes(x = x)) +
      geom_histogram(
        bins = 50,          
        fill = "steelblue",
        color = "white"
      ) +
      labs(
        title = bquote("Distribution of " * lambda[l] * " of " * .(name) * " in the " * .(celltype) * " tissue"),
        x = "Value",
        y = "Count"
      ) +
      theme_minimal()
    print(summary(limited))
    print(p)
    
  }
}


df_small <- cbind(lambda_result_small, Type = "Negative (< -1)")
df_large <- cbind(lambda_result_large, Type = "Positive (> 1)")

df_small$tissue <- rownames(df_small)
df_large$tissue <- rownames(df_large)
rownames(df_small) <- NULL
rownames(df_large) <- NULL
df_all <- rbind(df_small, df_large)
df_long <- df_all %>%
  pivot_longer(cols = c("Promoter", "Enhancer", "TFBS"),
               names_to = "Annotation",
               values_to = "Count") %>%
  rename(Tissue = tissue)

df_long$Tissue <- factor(df_long$Tissue,
                        levels = rev(rownames(lambda_result_small)))
df_long$Type <- factor(df_long$Type, levels = c("Negative (< -1)", "Positive (> 1)"))
df_long$Annotation <- factor(df_long$Annotation, levels = anno_name)


p = ggplot(df_long, aes(x = Count, y = Tissue, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Annotation, scales = "free_x") +
  scale_fill_manual(
    values = c("Negative (< -1)" = "#00BFC4",
               "Positive (> 1)" = "#F8766D"),
    labels = c(expression(Negative~(lambda[l] < -1)),
               expression(Positive~(lambda[l] > 1)))
  ) +
  labs(title = NULL,
       x = "Number of Genes", y = "Tissue") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 12, face = "bold"),
        legend.position = "top")


ggsave("/gpfs/gibbs/pi/zhao/xz527/annoKn_plots/lambda_l_type1.pdf", 
       p, 
       width = 10, height = 6, units = "in", 
       bg = "white", device = cairo_pdf)



################## check with the list of GTEx v10 eGenes: ##################

common1_all <- numeric()
common2_all <- numeric()
prop1_all <- numeric()
prop2_all <- numeric()


q_threshold = 0.05

allcelltype = c('Muscle_Skeletal', 'Whole_Blood',"Skin_Sun_Exposed_Lower_leg","Artery_Tibial",
                'Adipose_Subcutaneous', "Thyroid", "Nerve_Tibial", "Skin_Not_Sun_Exposed_Suprapubic")

for(celltype in allcelltype){
  
  print(celltype)
  egenes <- read.table(gzfile(paste0("~/GTEx/GTEx_Analysis_v10_eQTL_updated/",celltype,".v10.eGenes.txt.gz")),
                       header = TRUE,     
                       sep = "\t",        
                       stringsAsFactors = FALSE)
  
  head(egenes)
  egenes <- egenes[egenes$biotype == 'protein_coding',]
  egenes <- egenes[egenes$qval <= q_threshold,]
  
  
  sign1_t <- significant1_all[significant1_all$celltype == celltype, ]
  sign2_t <- significant2_all[significant2_all$celltype == celltype, ]
  
  common1 <- intersect(egenes$gene_name, sign1_t$hgnc_symbol)
  common2 <- intersect(egenes$gene_name, sign2_t$hgnc_symbol)
  
  common1_all <- append(common1_all, length(common1))
  common2_all <- append(common2_all, length(common2))
  
  print(length(common1))
  print(length(common2))
  
  prop1 <- length(common1) / nrow(sign1_t)
  prop2 <- length(common2) / nrow(sign2_t)
  
  prop1_all <- append(prop1_all, prop1)
  prop2_all <- append(prop2_all, prop2)
  
  cat("Proportion in sign1_t:", round(prop1, 3), "\n")
  cat("Proportion in sign2_t:", round(prop2, 3), "\n")
  
  
}

rbind(common1_all, common2_all)
rbind(prop1_all, prop2_all)




