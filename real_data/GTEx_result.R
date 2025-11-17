############ check the total number of risk genes:


pc = 5
# chrid = 1
lnc = FALSE
celltype = "Brain_Cortex"
number_gene <- numeric()

for(chrid in as.character(1:22)){
  print(chrid)
  library(data.table)
  library(genio)
  library(tidyr)
  library(dplyr)
  library(bigstatsr)
  
  
  bedNA <- function(bed1){
    for(j in 1:ncol(bed1)){
      temp <- bed1[,j]
      temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
      bed1[,j] <- temp
      #print(j)
    }
    return(bed1)
  }
  

  
  gexp_path = paste0('/gpfs/gibbs/pi/zhao/jz874/jiazhao/Xiangyu/TWAS_fm/GTEx_Gene_expression_adjust_covariates/adjusted_expr_age_sex_',pc,'genopcs_5peers')
  
  directory_path <- paste(gexp_path, "/chr",chrid, sep="")
  # List all folders in the specified directory
  folder_names <- list.dirs(directory_path, full.names = FALSE, recursive = FALSE)
  genes <- sub("\\..*", "", folder_names)
  
  # Print the folder names
  #print(folder_names)
  
  library(biomaRt)
  
  while (TRUE) {
    ensembl <- try({
      # useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
      useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
    })
    if (!isa(ensembl, "try-error"))
      break
    else{
      Sys.sleep(3)
      print(ensembl)
      print("try ensembl again!")
    }
  }
  
  print("useMart")
  
  while (TRUE) {
    z <- try({
      getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "transcript_biotype"),
            "ensembl_gene_id", genes, mart = ensembl)
    })
    if (!isa(z, "try-error"))
      break
    else{
      Sys.sleep(3)
      print(z)
      print("try z again!")
    }
    
  }
  
  index <- match(z$ensembl_gene_id, genes)
  z$folder_names <- folder_names[index]
  risk_gene <- z
  
  if(lnc == FALSE){
    risk_gene <- risk_gene[risk_gene$transcript_biotype %in% c("protein_coding"),]
  }else{
    risk_gene <- risk_gene[risk_gene$transcript_biotype %in% c("protein_coding", "lncRNA"),]
  }
  
  risk_gene <- risk_gene[!duplicated(risk_gene$ensembl_gene_id),]
  
  
  remain <- numeric()
  for(j in 1:nrow(risk_gene)){
    gene_info <- risk_gene[j,]
    path <- paste(gexp_path, "/chr",
                  chrid,"/",gene_info$folder_names,"/",celltype,".adj_expr", sep="")
    if (!file.exists(path)) {
      # cat("File does not exist for iteration", j, "\n")
      next  # Skip this iteration and move to the next one
    }
    remain <- append(remain, j)
  }
  
  risk_gene <- risk_gene[remain,]
  print(dim(risk_gene))
  number_gene <- append(number_gene, nrow(risk_gene))
}







celltype = "Whole_Blood"
celltype = "Muscle_Skeletal"
celltype = "Adipose_Subcutaneous"
celltype = "Artery_Tibial"
celltype = "Adipose_Visceral_Omentum"
celltype = "Brain_Cortex"
celltype = "Bladder"
celltype = "Pancreas"
celltype = "Liver"
celltype = "Thyroid"
celltype = "Skin_Sun_Exposed_Lower_leg"
celltype = "Nerve_Tibial"
celltype = "Breast_Mammary_Tissue"

knockoff.threshold <- function (W, fdr = 0.1, offset = 1){
  if (offset != 1 && offset != 0) {
    stop("Input offset must be either 0 or 1")
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t) (offset + sum(W <= -t))/max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}


all_tissue <- c("Whole_Blood", "Muscle_Skeletal", "Adipose_Subcutaneous", "Artery_Tibial",
                "Adipose_Visceral_Omentum", "Brain_Cortex", "Bladder", 
                "Liver", "Skin_Sun_Exposed_Lower_leg")

gene_all_tissue <- numeric()
knockoff_all_tissue <- numeric()
annoKn_all_tissue <- numeric()

for(celltype in all_tissue){
  
  
  gene_number <- c()
  
  sum1_all = c()
  sum2_all = c()
  
  
  find1 <- c()
  find2 <- c()
  lambdas <- c()
  
  
  chrid = '1'
  gene_i = '12'
  
  for(chrid in as.character(1:22)){
    
    alpha = 0.2
    
    sum1 = c()
    sum2 = c()
    
    for(gene_i in 1:2000){
      print(gene_i)
      # path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/chr_", chrid,"_gene_", gene_i, "_ct_", celltype, ".RData")
      path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/chr_", chrid,"_gene_", gene_i, "_ct_", celltype, ".RData")
      
      if (!file.exists(path)) {
        message("not exist: ", path)
        gene_number <- append(gene_number, gene_i - 1)
        break
      }
      
      load(path)
      
      W1 = final_result$ghostknockoff
      tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
      rej1 = as.numeric(which(W1>=tau))
      # print(rej1)
      
      
      result_raw <- final_result$knockoff_anno
      
      beta = result_raw$beta
      p = length(beta)/2
      T = abs(beta[1:p])
      T_tilde = abs(beta[(p+1):(2*p)])
      T_max = pmax(T,T_tilde)
      W = T - T_tilde
      
      # hist(W, breaks = 50)
      tau = knockoff.threshold(W,fdr = alpha,offset = 1)
      rej2 = as.numeric(which(W >= tau))
      
      result_raw$lambda_s
      print("rej1:")
      print(rej1)
      if(length(rej1 > 0)){
        sum1 = append(sum1, gene_i)
        find1 <- append(find1, paste(chrid, gene_i))
      }
      print("rej2:")
      print(rej2)
      if(length(rej2 > 0)){
        sum2 = append(sum2, gene_i)
        find2 <- append(find2, paste(chrid, gene_i))
        lambdas <- rbind(lambdas, final_result$knockoff_anno$lambda_s[[1]])
      }
      
    }
    sum1_all <- append(sum1_all, length(sum1))
    sum2_all <- append(sum2_all, length(sum2))
    
  }
  
  gene_all_tissue <- append(gene_all_tissue, sum(gene_number))
  
  knockoff_all_tissue <- append(knockoff_all_tissue, length(find1))
  annoKn_all_tissue <- append(annoKn_all_tissue, length(find2))
  
  
  library(grid)
  grid.newpage()
  library(VennDiagram)
  
  venn.plot <- draw.pairwise.venn(
    area1 = length(find1),                 
    area2 = length(find2),                 
    cross.area = length(intersect(find1, find2)),  
    category = c("Knockoff", "AnnoKn"),
    fill = c("skyblue", "orange"),
    lty = "blank",
    cex = 1.5,
    cat.cex = 1.5
  )
  
  grid.draw(venn.plot)
  
}


data.frame(tissue = all_tissue, 
           knockoff = knockoff_all_tissue, 
           annoKn = annoKn_all_tissue,
           improve = (annoKn_all_tissue - knockoff_all_tissue)/knockoff_all_tissue)




labels <- paste0("chr", 1:22)

df <- data.frame(
  sum1 = sum1_all,
  sum2 = sum2_all,
  chr  = labels
)

ggplot(df, aes(x = sum1, y = sum2, label = chr)) +
  geom_point(color = "blue", size = 3) +
  geom_text(vjust = -0.5, color = "red") +
  xlim(c(0, 160)) +
  ylim(c(0, 160)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Knockoff", y = "Knockoff-anno",
       title = "Scatter Plot with Chromosome Labels") +
  theme_minimal()












sum_all <- numeric()
sum4_all = c()
sum3_all = c()
find4 <- c()
find3 <- c()
for(chrid in 1:22){
  sum4 = c()
  sum3 = c()
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/adakn_chr_", chrid, "_ct_", celltype, ".RData")
  
  if (!file.exists(path)) {
    message("not exist: ", path)
    break
  }
  load(path)
  
  print(length(all_result))
  
  sum_all <- append(sum_all, length(all_result))
  
  for(gene_i in 1:length(all_result)){
    print(gene_i)
    adakn_result = all_result[[gene_i]]
    
    W1 = adakn_result$ghostknockoff
    tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
    rej1 = as.numeric(which(W1>=tau))
    print("rej1:")
    print(rej1)
    
    if(length(rej1 > 0)){
      sum4 = append(sum4, gene_i)
      find4 <- append(find4, paste(chrid, gene_i))
    }
    sum4_all <- append(sum4_all, length(sum4))
    
    if(length(adakn_result$adakn) == 1){
      rej3 = rej1
    }else{
      rej3 = adakn_result$adakn$rejs[[2]]
    }
    
    print("rej3:")
    print(rej3)
    if(length(rej3 > 0)){
      sum3 = append(sum3, gene_i)
      find3 <- append(find3, paste(adakn_result$chr, adakn_result$gene_i))
    }
    sum3_all <- append(sum3_all, length(sum3))
  }
}





celltype = 'Whole_Blood'
celltype = "Brain_Cortex"
celltype = 'Muscle_Skeletal'
celltype = 'Adipose_Subcutaneous'
celltype = "Artery_Tibial"
celltype = "Adipose_Visceral_Omentum"
celltype = "Lung"
celltype = 'Pancreas'
celltype = 'Liver'
celltype = "Skin_Sun_Exposed_Lower_leg"
celltype = "Thyroid"
celltype = "Nerve_Tibial"
celltype = "Breast_Mammary_Tissue"



celltype = "Brain_Cortex"



########### new result

sum_all_ct <- data.frame(matrix(nrow = 22, ncol = 0))
number_gene_ct <- data.frame(matrix(nrow = 22, ncol = 0))
find1_all <- list()
find2_all <- list()



allcelltype = c('Muscle_Skeletal', 'Whole_Blood',"Skin_Sun_Exposed_Lower_leg","Artery_Tibial",
                'Adipose_Subcutaneous', "Thyroid", "Nerve_Tibial", "Skin_Not_Sun_Exposed_Suprapubic", "Lung", "Esophagus_Mucosa")

#, "Thyroid", "Nerve_Tibial", "Skin_Not_Sun_Exposed_Suprapubic",
# "Lung", "Esophagus_Mucosa"

type = 3
alpha = 0.2
pc = 5
lnc = FALSE



library(data.table)
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


for(celltype in allcelltype){
  
  lambda_s <- numeric()
  
  number_gene <- numeric()
  
  for(chrid in as.character(1:22)){
    print(celltype)
    print(chrid)
    library(data.table)
    library(genio)
    library(tidyr)
    library(dplyr)
    library(bigstatsr)
    
    
    bedNA <- function(bed1){
      for(j in 1:ncol(bed1)){
        temp <- bed1[,j]
        temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
        bed1[,j] <- temp
        #print(j)
      }
      return(bed1)
    }
    
    
    
    gexp_path = paste0('/gpfs/gibbs/pi/zhao/jz874/jiazhao/Xiangyu/TWAS_fm/GTEx_Gene_expression_adjust_covariates/adjusted_expr_age_sex_',pc,'genopcs_5peers')
    
    directory_path <- paste(gexp_path, "/chr",chrid, sep="")
    # List all folders in the specified directory
    folder_names <- list.dirs(directory_path, full.names = FALSE, recursive = FALSE)
    genes <- sub("\\..*", "", folder_names)
    
    library(biomaRt)
    
    while (TRUE) {
      ensembl <- try({
        # useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
        useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
      })
      if (!isa(ensembl, "try-error"))
        break
      else{
        Sys.sleep(3)
        print(ensembl)
        print("try ensembl again!")
      }
    }
    
    print("useMart")
    
    while (TRUE) {
      z <- try({
        getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "transcript_biotype"),
              "ensembl_gene_id", genes, mart = ensembl)
      })
      if (!isa(z, "try-error"))
        break
      else{
        Sys.sleep(3)
        print(z)
        print("try z again!")
      }
      
    }
    
    index <- match(z$ensembl_gene_id, genes)
    z$folder_names <- folder_names[index]
    risk_gene <- z
    
    if(lnc == FALSE){
      risk_gene <- risk_gene[risk_gene$transcript_biotype %in% c("protein_coding"),]
    }else{
      risk_gene <- risk_gene[risk_gene$transcript_biotype %in% c("protein_coding", "lncRNA"),]
    }
    
    risk_gene <- risk_gene[!duplicated(risk_gene$ensembl_gene_id),]
    
    
    remain <- numeric()
    for(j in 1:nrow(risk_gene)){
      gene_info <- risk_gene[j,]
      path <- paste(gexp_path, "/chr",
                    chrid,"/",gene_info$folder_names,"/",celltype,".adj_expr", sep="")
      if (!file.exists(path)) {
        next  # Skip this iteration and move to the next one
      }
      remain <- append(remain, j)
    }
    
    risk_gene <- risk_gene[remain,]
    print(dim(risk_gene))
    number_gene <- append(number_gene, nrow(risk_gene))
  }
  
  
  
  sum_all <- numeric()
  sum1_all = c()
  sum2_all = c()
  find1 <- c()
  find2 <- c()
  for(chrid in 1:22){
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
      
      W1 = adakn_result$ghostknockoff
      tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
      rej1 = as.numeric(which(W1>=tau))
      # print("rej1:")
      # print(rej1)
      if(length(rej1 > 0)){
        sum1 = append(sum1, gene_i)
        find1 <- append(find1, paste(chrid, gene_i))
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
      }
    }
    sum1_all <- append(sum1_all, length(sum1))
    sum2_all <- append(sum2_all, length(sum2))
  }
  
  cbind(number_gene, sum_all)
  
  print(length(find1))
  print(length(find2))
  
  sum_all_ct[[celltype]] = sum_all
  number_gene_ct[[celltype]] = number_gene
  
  find1_all <- append(find1_all, list(find1))
  find2_all <- append(find2_all, list(find2))
  
  lambdas_all <- append(lambdas_all, list(lambda_s))
  
}

sum_all_ct_0 = sum_all_ct
number_gene_ct_0 = number_gene_ct

number_gene_ct - sum_all_ct
# a = numeric()
# for(adakn_result in all_result){
#   a = append(a, adakn_result$gene_i)
# }


improve_prop_0 = (sapply(find2_all, length) - sapply(find1_all, length))/sapply(find1_all, length)



################################################ lambda_s

i = 4
j = 3

celltype = allcelltype[i]
name = c("Promoter", "Enhancer", "TFBS")[j]

limited <- lambdas_all[[i]][, j]
limited <- limited[which(abs(limited) < 2)]
df <- data.frame(x = limited)

ggplot(df, aes(x = x)) +
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

summary(limited)

















################## check the total number of protein coding genes


number_gene_ct <- data.frame(matrix(nrow = 22, ncol = 0))


allcelltype = c('Muscle_Skeletal', 'Whole_Blood',"Skin_Sun_Exposed_Lower_leg","Artery_Tibial",
                'Adipose_Subcutaneous', "Thyroid", "Nerve_Tibial", "Skin_Not_Sun_Exposed_Suprapubic", "Lung", "Esophagus_Mucosa")

type = 1
alpha = 0.1
pc = 5
lnc = FALSE

lambdas_all <- list()


for(celltype in allcelltype){
  
  lambda_s <- numeric()
  
  number_gene <- numeric()
  
  for(chrid in as.character(1:22)){
    print(celltype)
    print(chrid)
    library(data.table)
    library(genio)
    library(tidyr)
    library(dplyr)
    library(bigstatsr)
    
    
    bedNA <- function(bed1){
      for(j in 1:ncol(bed1)){
        temp <- bed1[,j]
        temp[is.na(temp)] <- mean(temp,na.rm = TRUE)
        bed1[,j] <- temp
        #print(j)
      }
      return(bed1)
    }
    
    
    
    gexp_path = paste0('/gpfs/gibbs/pi/zhao/jz874/jiazhao/Xiangyu/TWAS_fm/GTEx_Gene_expression_adjust_covariates/adjusted_expr_age_sex_',pc,'genopcs_5peers')
    
    directory_path <- paste(gexp_path, "/chr",chrid, sep="")
    # List all folders in the specified directory
    folder_names <- list.dirs(directory_path, full.names = FALSE, recursive = FALSE)
    genes <- sub("\\..*", "", folder_names)
    
    library(biomaRt)
    
    while (TRUE) {
      ensembl <- try({
        # useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
        useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host="https://www.ensembl.org")
      })
      if (!isa(ensembl, "try-error"))
        break
      else{
        Sys.sleep(3)
        print(ensembl)
        print("try ensembl again!")
      }
    }
    
    print("useMart")
    
    while (TRUE) {
      z <- try({
        getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "transcript_biotype"),
              "ensembl_gene_id", genes, mart = ensembl)
      })
      if (!isa(z, "try-error"))
        break
      else{
        Sys.sleep(3)
        print(z)
        print("try z again!")
      }
      
    }
    
    index <- match(z$ensembl_gene_id, genes)
    z$folder_names <- folder_names[index]
    risk_gene <- z
    
    if(lnc == FALSE){
      risk_gene <- risk_gene[risk_gene$transcript_biotype %in% c("protein_coding"),]
    }else{
      risk_gene <- risk_gene[risk_gene$transcript_biotype %in% c("protein_coding", "lncRNA"),]
    }
    
    risk_gene <- risk_gene[!duplicated(risk_gene$ensembl_gene_id),]
    
    
    remain <- numeric()
    for(j in 1:nrow(risk_gene)){
      gene_info <- risk_gene[j,]
      path <- paste(gexp_path, "/chr",
                    chrid,"/",gene_info$folder_names,"/",celltype,".adj_expr", sep="")
      if (!file.exists(path)) {
        next  # Skip this iteration and move to the next one
      }
      remain <- append(remain, j)
    }
    
    risk_gene <- risk_gene[remain,]
    print(dim(risk_gene))
    number_gene <- append(number_gene, nrow(risk_gene))
  }
  
  number_gene_ct[[celltype]] = number_gene
}

