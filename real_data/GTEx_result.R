############ check the total number of risk genes:


pc = 5
# chrid = 1
lnc = FALSE
celltype = "Whole_Blood"
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


gene_number <- c()
sum1_all = c()
sum2_all = c()

find1 <- c()
find2 <- c()

chrid = '1'
gene_i = '12'

for(chrid in as.character(1:22)){
  
  alpha = 0.2
  
  sum1 = c()
  sum2 = c()
  
  
  for(gene_i in 1:2000){
    print(gene_i)
    path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/GTEx/result_new_500k_enhancer/chr_", chrid,"_gene_", gene_i, "_ct_", celltype, ".RData")
    
    if (!file.exists(path)) {
      message("not exist: ", path)
      break
    }
    gene_number <- append(gene_number, gene_i - 1)
    
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
    }
  }
  sum1_all <- append(sum1_all, length(sum1))
  sum2_all <- append(sum2_all, length(sum2))
}


library(grid)
grid.newpage()
library(VennDiagram)

venn.plot <- draw.pairwise.venn(
  area1 = length(find1),                 
  area2 = length(find2),                 
  cross.area = length(intersect(find1, find2)),  
  category = c("Knockoff", "Knockoff-anno"),
  fill = c("skyblue", "orange"),
  lty = "blank",
  cex = 1.5,
  cat.cex = 1.5
)

grid.draw(venn.plot)










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
