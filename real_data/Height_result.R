
########## check the risk region for different CHRs

s = numeric()

for(chrid in 1:22){
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/Height/annot_pvalus/risk_region/","EUR","_chr_", chrid,".txt")
  risk_region <- read.table(path)
  if(any(is.na(risk_region))){
    s = append(s, 0)
    next
  }
  print(dim(risk_region))
  s = append(s, nrow(risk_region))
}


chr <- 1:length(s)
df <- data.frame(Chromosome = factor(chr), Count = s)

library(ggplot2)

ggplot(df, aes(x = Chromosome, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Number of Risk Regions per Chromosome",
    x = "Chromosome",
    y = "Risk Region Count"
  )

########################

region_1 <- numeric()
region_anno_1 <- numeric()
region_2 <- numeric()
region_anno_2 <- numeric()
lambdas1 <- numeric()
lambdas2 <- numeric()

threshold = 0.1
sum1 = 0
sum2 = 0
chrid = '5'
M = 1
seed = '1000'
# seed = '12345' for dss


chrlist = (1:22)
for(chrid in chrlist){
  print(paste0("chrid", chrid))
  ancestry0 = 'EUR'
  ancestry1 = c('AFR')
  other = paste(ancestry1, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/Height/annot_pvalus/result/", ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)
  
  print(length(result$risk_region))
  
  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)
  
  
  
  for(i in 1:n_region){
    # print(i)
    a = q_gk[[i]]
    if(sum(a <= threshold) > 0){
      region_1 <- rbind(region_1, as.numeric(append(chrid, result$risk_region[[i]])))
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      region_anno_1 <- rbind(region_anno_1, as.numeric(append(chrid, result$risk_region[[i]])))
      lambdas1 <- append(lambdas1, result$lambda_s[[i]])
    }
    
  }
  
  
  ancestry0 = 'AFR'
  ancestry1 = c('EUR')
  other = paste(ancestry1, collapse = "_")
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/real_data/Height/annot_pvalus/result/", ancestry0, "_result_final_10_median_r2_", other, "_chr_",chrid, "_M_", M, "_", seed, ".RData")
  load(path)
  print(length(result$risk_region))
  
  q_gk = result$q_gk
  q_gkanno = result$q_gkanno
  n_region = length(q_gk)
  
  
  for(i in 1:n_region){
    # print(i)
    a = q_gk[[i]]
    if(sum(a <= threshold) > 0){
      region_2 <- rbind(region_2, as.numeric(append(chrid, result$risk_region[[i]])))
    }
    b = q_gkanno[[i]]
    # print(which(b <= threshold))
    if(sum(b <= threshold) > 0){
      region_anno_2 <- rbind(region_anno_2, as.numeric(append(chrid, result$risk_region[[i]])))
      lambdas2 <- append(lambdas2, result$lambda_s[[i]])
    }
    
  }
  
}

nrow(region_1)
nrow(region_2)

nrow(region_anno_1)
nrow(region_anno_2)




df1 <- as.data.frame(region_1)
df2 <- as.data.frame(region_2)
common_rows <- merge(df1, df2)
n_common <- nrow(common_rows)
n_common


df1 <- as.data.frame(region_anno_1)
df2 <- as.data.frame(region_anno_2)
common_rows <- merge(df1, df2)
n_common <- nrow(common_rows)
n_common




library(ggplot2)
df <- data.frame(lambdas = lambdas1)
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,  
    fill = "#0072B2", 
    color = "white",
    alpha = 0.8     
  ) +
  labs(
    title = expression("Distribution of " * lambda[l] * " of AFR p-values for EUR GWAS"),
    x = expression(lambda[l]),         
    y = "Frequency (Count)"         
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot



df <- data.frame(lambdas = lambdas2)
histogram_plot <- ggplot(df, aes(x = lambdas)) +
  geom_histogram(
    bins = 30,  
    fill = "#0072B2", 
    color = "white",
    alpha = 0.8     
  ) +
  labs(
    title = expression("Distribution of " * lambda[l] * " of EUR p-values for AFR GWAS"),
    x = expression(lambda[l]),       
    y = "Frequency (Count)"         
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 12)
  )
histogram_plot




library(VennDiagram)
library(grid)

df1 <- as.data.frame(region_2)
df2 <- as.data.frame(region_anno_2)
set1 <- apply(df1, 1, paste, collapse = "_")
set2 <- apply(df2, 1, paste, collapse = "_")

n1 <- length(set1)
n2 <- length(set2)
n_common <- length(intersect(set1, set2))

venn.plot <- draw.pairwise.venn(
  area1 = n1,
  area2 = n2,
  cross.area = n_common,
  category = c("GhostKnockoff", "AnnoGK"),
  fill = c("steelblue", "orange"),
  
  alpha = c(0.6, 0.6),
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

grid.newpage()
grid.draw(venn.plot)

grid.text(
  "Risk regions identified by GhostKnockoff and AnnoGK in African",
  x = 0.5, y = 0.98,          
  gp = gpar(fontsize = 16, fontface = "bold")  
)

