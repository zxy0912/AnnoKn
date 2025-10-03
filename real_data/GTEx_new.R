
args = commandArgs(trailingOnly=TRUE)
options(stringsAsFactors=F)

chrid = as.numeric(args[1])

source("/home/xz527/Rcode/knockoff_anno/KF_anno/KF_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GK_anno.R")
source("/home/xz527/Rcode/knockoff_anno/GK_anno/GhostKnockoff.R")
packageVersion("Matrix") 

################################## package for AdaKn ##########################################
suppressMessages(library("adaptMT"))
suppressMessages(library("splines"))
suppressMessages(library("knockoff"))
suppressMessages(library("SNPknock"))
suppressMessages(library("dplyr"))
suppressMessages(library("corpcor"))
suppressMessages(library("glmnet"))
suppressMessages(library("MASS"))
suppressMessages(library("gam"))
suppressMessages(library("randomForest"))
suppressMessages(library("tidyverse"))
suppressMessages(library("mgcv"))
source_gitfile <- function(filename){
  source(sprintf("%s.R", filename))
}
file_vec <- c("/home/xz527/Rcode/knockoff_side/code//utils/all_other_methods",
              "/home/xz527/Rcode/knockoff_side/code/utils/adaptive_knockoff",
              "/home/xz527/Rcode/knockoff_side/code/utils/filter_EM",
              "/home/xz527/Rcode/knockoff_side/code/utils/filter_gam",
              "/home/xz527/Rcode/knockoff_side/code/utils/filter_glm",
              "/home/xz527/Rcode/knockoff_side/code/utils/filter_randomForest", 
              "/home/xz527/Rcode/knockoff_side/code/utils/All_q_est_functions", 
              "/home/xz527/Rcode/knockoff_side/code/utils/accumulation_test_functions")
getfile <- sapply(file_vec,source_gitfile)
############################################################################


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


pc = 5
# chrid = 1
lnc = FALSE
celltype = "Muscle_Skeletal"
# "Brain_Cortex" "Lung" "Pancreas"

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
    cat("File does not exist for iteration", j, "\n")
    next  # Skip this iteration and move to the next one
  }
  remain <- append(remain, j)
}

risk_gene <- risk_gene[remain,]



########### read in the genotype data and match to snp

maf = 0.01

if(maf == 0.05){
  maf_indi = 'maf5'
  maf_dir = 'bychr_maf0.05'
}else if(maf == 0.01){
  maf_indi = 'maf1'
  maf_dir = 'bychr_maf0.01'
}

ancestry = 'EUR'
path <- paste("/gpfs/gibbs/pi/zhao/xz527/TWAS_fm/real_data/GTEX/ancestry/",ancestry,"/", maf_dir, "/GTEX_chr", chrid,sep='')
gtex_chr <- read_plink(path)

snp_mart <- useEnsembl(
  biomart = "ENSEMBL_MART_SNP",
  dataset = "hsapiens_snp"
)



cats = c("SNP", "Promoter_UCSC", "Enhancer_Hoffman", "Coding_UCSC", "TFBS_ENCODE")
cats = c("SNP", "Promoter_UCSC", "Enhancer_Hoffman", "TFBS_ENCODE")
# cats = c("SNP","Coding_UCSC", "TSS_Hoffman", "Enhancer_Hoffman", 
#          "Promoter_UCSC", "H3K27ac_Hnisz", "H3K4me1_peaks_Trynka")
# annot = fread(paste0("/gpfs/gibbs/pi/zhao/cl2384/MESC_Genes/baselineLD_2.1/baselineLD.", chrid, ".annot.gz"))
annot_all = fread(paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/LDSC/data/GRCh38/baselineLD_v2.2/baselineLD.", chrid, ".annot.gz"))
annot_all = annot_all[, ..cats]

######################################################################################
# path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/ucsc_hg38/chr",chrid,".txt")
# ref_hg38 <- read.table(path, header = FALSE)
# colnames(ref_hg38) <- c("chr", "start","end","snp")
# 
# ref_hg38_selected <- ref_hg38[ref_hg38$end %in% gtex_chr$bim$pos | ref_hg38$snp %in% annot_all$SNP,]
# path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/ucsc_hg38/chr",chrid,"_selected.txt")
# write.table(ref_hg38_selected, file = path, row.names = FALSE,
#             col.names = TRUE, quote = FALSE)
######################################################################################

path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/ucsc_hg38/chr",chrid,"_selected.txt")
ref_hg38 <- read.table(path, header = TRUE)


####################
# filtered_ref <- ref_hg38 %>%
#   group_by(end) %>%
#   slice(if (any(snp %in% annot_all$SNP)) {
#     which(snp %in% annot_all$SNP)[1]
#   } else {
#     sample(1:n(), 1)
#   }) %>%
#   ungroup()

##############

M = 250000

gene_i = 1

for(gene_i in 1:nrow(risk_gene)){
  
  set.seed(1000)
  
  gene_info <- risk_gene[gene_i,]
  path <- paste(gexp_path, "/chr",
                chrid,"/",gene_info$folder_names,"/",celltype,".adj_expr", sep="")
  
  start = gene_info$start_position - M
  end = gene_info$end_position + M
  
  # check
  if(j==1){
    print(paste("read gene expression data from", path))
  }
  
  expr <- fread(path)
  expr = as.data.frame(expr)
  short_ids <- sapply(strsplit(expr$individual, "-"), function(x) paste(x[1:2], collapse = "-"))
  expr$individual = short_ids
  colnames(expr) = c("V1","V2")
  
  
  
  ######### snps that within this gene:
  index <- which(gtex_chr$bim$pos > start & gtex_chr$bim$pos < end)
  # dim(gtex_chr$X)
  X_gtex_region <- gtex_chr$X[index, ]
  
  individuals <- colnames(X_gtex_region)
  common_ind <- intersect(expr$V1, individuals)
  gidx <- match(common_ind, expr$V1)
  midx <- match(common_ind, individuals)
  
  expr <- expr[gidx,]
  y = expr$V2
  
  if(var(y) == 0){
    cat("variance of gene expression is 0:",
        "  chrid:", chrid, 
        "  gene:", gene_i, "\n", 
        file = "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/log_", celltype, ".txt", append = TRUE)
    next
  }
  
  
  ####### remove constant SNPs 
  X_region <- X_gtex_region[,midx]
  maf <- apply(X_region, 1, function(x){mean(x, na.rm = TRUE)})
  bad_snp = as.numeric(which(apply(X_region, 1, var)==0))
  
  if(length(bad_snp) > 0){
    print(paste0("gene:", j, "  bad_snp:", bad_snp))
    # snpindex = snpindex[-bad_snp]
    X_region <- X_region[-bad_snp,]
  }
  
  ######## cluster and select representative SNPs
  
  X_region <- t(bedNA(X_region))
  X_region <- scale(X_region)
  
  table(rownames(X_region) == expr$V1)
  
  ########## calculate the pvalue
  
  X_fbm <- as_FBM(X_region)
  res <- big_univLinReg(X_fbm, y.train = y)
  pvals <- predict(res, log10 = FALSE)
  
  ################
  
  
  R = cor(X_region)
  d <- as.dist(1 - R)  # higher distance = lower correlation
  hc <- hclust(d, method = "average")
  plot(hc, labels = FALSE, main = "SNP clustering by LD")
  clusters <- cutree(hc, h = 0.5)  # e.g., cluster SNPs with r² > 0.5
  length(table(clusters))
  mean(table(clusters)) # average size for each cluster
  
  if(length(table(clusters)) < 11){
    cat("number of SNP groups <= 10:",
        "  chrid:", chrid, 
        "  gene:", gene_i, "\n", 
        file = "/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/log_", celltype, ".txt", append = TRUE)
    next
  }
  
  
  df <- data.frame(
    SNP = names(clusters),
    pval = pvals,
    cluster = clusters,
    index = 1:length(pvals)
  )
  
  representatives <- df %>%
    group_by(cluster) %>%
    slice_min(order_by = pval, n = 1, with_ties = FALSE)
  
  head(representatives)
  dim(representatives)
  
  
  index <- representatives$index[order(representatives$index)]
  
  X_repre = X_region[, index]
  y_repre = (y - mean(y))/sd(y)
  
  
  print(dim(X_repre))
  
  
  
  #################### convert to rsid
  
  name_snp <- colnames(X_repre)
  df <- data.frame(name_snp = name_snp)
  df_split <- df %>%
    separate(name_snp, into = c("chr", "pos", "ref", "alt", "build"), sep = "_") %>%
    mutate(chr = gsub("^chr", "", chr),  
           pos = as.integer(pos))
  
  head(df_split)
  # position <- paste(df_split$chr, df_split$pos, df_split$pos, sep = ":")
  
  
  # while (TRUE) {
  #   snp_info <- try({
  #     getBM(
  #       attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele'),
  #       filters = 'chromosomal_region',
  #       values = position,
  #       mart = snp_mart
  #     )
  #   })
  #   if (!isa(snp_info, "try-error"))
  #     break
  #   else{
  #     break
  #     Sys.sleep(3)
  #     # print(z)
  #     print("try again!")
  #   }
  # }
  # 
  # if (!isa(snp_info, "try-error")){
  #   print("try error")
  #   next
  # }
  
  # index <- match(df_split$pos, snp_info$chrom_start)
  # snp_name <- snp_info$refsnp_id[index]
  # colnames(X_repre) <- snp_name
  
  
  
  length(intersect(df_split$pos, ref_hg38$end))
  index <- match(df_split$pos, ref_hg38$end)
  snp_name <- ref_hg38$snp[index]
  
  ############################################################
  ref_hg38_temp <- ref_hg38[ref_hg38$end %in% df_split$pos,]
  ref_hg38_temp <- ref_hg38_temp %>%
    group_by(end) %>%
    slice(if (any(snp %in% annot_all$SNP)) {
      which(snp %in% annot_all$SNP)[1]
    } else {
      sample(1:n(), 1)
    }) %>%
    ungroup()
  
  index <- match(df_split$pos, ref_hg38_temp$end)
  snp_name <- ref_hg38_temp$snp[index]
  ############################################################
  
  
  # is snp not exists in ref_hg38, which is very unlikely, name this SNP the position on this chromosome
  snp_name[is.na(snp_name)] <- df_split$pos[is.na(snp_name)]
  colnames(X_repre) <- snp_name
  
  
  
  
  
  
  ############## knockoff
  
  # set.seed(111)
  Xk = create.second_order(X_repre, method = 'sdp')
  
  # standardize the matrix
  Xk = scale(Xk)
  
  # run original knockoff
  p = ncol(X_repre)
  mdl = cv.glmnet(cbind(X_repre, Xk), y_repre, alpha=1)
  cvlambda = mdl$lambda.min
  beta = mdl$glmnet.fit$beta[, mdl$lambda == cvlambda]
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W1 = T-T_tilde
  as.numeric(W1)
  hist(W1, breaks = 50)
  
  alpha = 0.2
  tau = knockoff.threshold(W1,fdr = alpha,offset = 1)
  rej1 = as.numeric(which(W1>=tau))
  print(rej1)
  
  
  
  ############## knockoff with annotation
  # "Enhancer_Andersson", "GTEx_eQTL_MaxCPP"
  
  length(intersect(colnames(X_repre), annot_all$SNP))
  index <- match(colnames(X_repre), annot_all$SNP)
  annot <- annot_all[index, ]
  annot0 <- annot
  annot <- annot[, -1]
  
  apply(annot, 2, table)
  
  R <- apply(annot, 2, function(col){
    col[is.na(col)] <- mean(col, na.rm = TRUE)
    col
  })
  
  
  ############################### AdaKn ################################################
  alphalist <- c(0.05, 0.1, 0.2, 0.3)
  res_result = filter_randomForest(W1,R,alpha =alphalist,offset=1)
  rej3 = res_result$rejs[[2]]
  rej3
  
  ######################################################################################
  
  # var_para <- apply(R, 2, var)
  # if(length(var_para!=0) == 0){
  #   print("no useful annotations")
  #   R = matrix(1, nrow = nrow(R), ncol = 1)
  # }else{
  #   R <- R[, var_para!=0]
  # }
  # 
  # R = scale(R)
  
  
  result_raw <- knockoff_anno_improved(X = X_repre, Xk = Xk, y = y_repre, attempts = c(0), R = R)
  
  beta = result_raw$beta
  T = abs(beta[1:p])
  T_tilde = abs(beta[(p+1):(2*p)])
  T_max = pmax(T,T_tilde)
  W = T - T_tilde
  
  hist(W, breaks = 50)
  
  
  tau = knockoff.threshold(W,fdr = alpha,offset = 1)
  rej2 = as.numeric(which(W >= tau))
  
  result_raw$lambda_s
  
  # hist(W[W!=0],  breaks = 20)
  # hist(W1[W1!=0],  breaks = 20)
  # hist(result_raw$Prior_final, xlab = 'weight', main = 'Histogram of weights')
  
  
  print(result_raw$lambda_s)
  print(rej1)
  print(rej2)
  
  
  final_result = list(ghostknockoff = W1,
                      knockoff_anno = result_raw,
                      adakn = res_result,
                      snp_name = snp_name,
                      annot = annot0,
                      R = R,
                      gene_info = gene_info)
  
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/adakn_chr_", chrid,"_gene_", gene_i, "_ct_", celltype, ".RData")
  # path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/chr_", chrid,"_gene_", gene_i, "_ct_", celltype, ".RData")
  save(final_result, file = path)
  
  # path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result/chr_", chrid,"_gene_", gene_i, ".RData")
  # cats = c("SNP", "Promoter_UCSC", "Enhancer_Hoffman", "Coding_UCSC", "TFBS_ENCODE")
  # M = 500000
  # clusters <- cutree(hc, h = 0.5) 
}

