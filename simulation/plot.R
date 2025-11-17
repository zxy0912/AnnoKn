
library(ggplot2)

iteration = 100

alphalist <- seq(0.3,0.05,-0.05)
len = length(alphalist)

############################ 1 dimension annotation ############################## 

power1 = matrix(0, nrow = iteration, ncol = len)
power2 = matrix(0, nrow = iteration, ncol = len)
power3 = matrix(0, nrow = iteration, ncol = len)
power4 = matrix(0, nrow = iteration, ncol = len)
power5 = matrix(0, nrow = iteration, ncol = len)
power6 = matrix(0, nrow = iteration, ncol = len)
power7 = matrix(0, nrow = iteration, ncol = len)
power8 = matrix(0, nrow = iteration, ncol = len)
fdr1 = matrix(0, nrow = iteration, ncol = len)
fdr2 = matrix(0, nrow = iteration, ncol = len)
fdr3 = matrix(0, nrow = iteration, ncol = len)
fdr4 = matrix(0, nrow = iteration, ncol = len)
fdr5 = matrix(0, nrow = iteration, ncol = len)
fdr6 = matrix(0, nrow = iteration, ncol = len)
fdr7 = matrix(0, nrow = iteration, ncol = len)
fdr8 = matrix(0, nrow = iteration, ncol = len)

s = 0
removelist = c()

binary = 'binary'
dimen = 'double'
amp = 3.5

for (i in 1:iteration) {
  path = paste0("/gpfs/gibbs/project/zhao/cl2667/knockoff_side/iter_simu/double/1dimS/result_", binary, "_amp_",amp, "_", i,".RData")
  #path = paste0("/gpfs/gibbs/project/zhao/cl2667/knockoff_side/iter_simu/double/1dimSo/result_", binary, "_amp_",amp, "_", i,".RData")
  if (!file.exists(path)){
    print(paste("iteration",i,"doesn't exist"))
    removelist = append(removelist, i)
    next
  }
  load(path) 
  s  = s+1
  
  power1[i,] <- result$power1
  power2[i,] <- result$power2
  power3[i,] <- result$power3
  power4[i,] <- result$power4
  power5[i,] <- result$power5
  power6[i,] <- result$power6
  power7[i,] <- result$power7
  power8[i,] <- result$power8
  
  
  fdr1[i,]   <- result$fdr1
  fdr2[i,]   <- result$fdr2
  fdr3[i,]   <- result$fdr3
  fdr4[i,]   <- result$fdr4
  fdr5[i,]   <- result$fdr5
  fdr6[i,]   <- result$fdr6
  fdr7[i,]   <- result$fdr7
  fdr8[i,]   <- result$fdr8
  
}

if(length(removelist) > 0){
  power1 = power1[-removelist,]
  power2 = power2[-removelist,]
  power3 = power3[-removelist,]
  power4 = power4[-removelist,]
  fdr1 = fdr1[-removelist,]
  fdr2 = fdr2[-removelist,]
  fdr3 = fdr3[-removelist,]
  fdr4 = fdr4[-removelist,]
}


power <- c(
  colMeans(power1, na.rm = TRUE),
  colMeans(power2, na.rm = TRUE),
  #colMeans(power3, na.rm = TRUE),
  #colMeans(power4, na.rm = TRUE),
  colMeans(power5, na.rm = TRUE),
  colMeans(power6, na.rm = TRUE),
  colMeans(power7, na.rm = TRUE),
  colMeans(power8, na.rm = TRUE)
)


fdr <- c(
  colMeans(fdr1, na.rm = TRUE),
  colMeans(fdr2, na.rm = TRUE),
  #colMeans(fdr3, na.rm = TRUE),
  #colMeans(fdr4, na.rm = TRUE),
  colMeans(fdr5, na.rm = TRUE),
  colMeans(fdr6, na.rm = TRUE),
  colMeans(fdr7, na.rm = TRUE),
  colMeans(fdr8, na.rm = TRUE)
)

method <- rep(c("AdaPT","Knockoff","AdaKn(RF)", "AdaKn(EM)", "AnnoKn", "AnnoKn-lite"), each = len)
plot_data <- data.frame(alpha = rep(alphalist, 6),
                        method = method,
                        power = power,
                        fdr = fdr)


p1 <- ggplot(plot_data, aes(x = alpha, y = power, color = method)) +
  geom_line() +
  geom_point() +
  labs(x = "Target FDR", y = "Power") +
  theme_minimal()

p2 <- ggplot(plot_data, aes(x = alpha, y = fdr, color = method)) +
  geom_line() +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Target FDR", y = "FDR") +
  theme_minimal()

cols <- c(
  "Knockoff"    = "#E69F00",
  "AnnoKn-lite" = "#56B4E9",
  "AnnoKn"      = "#FF0000", 
  "AdaKn(RF)"   = "#D55E00",
  "AdaKn(EM)"   = "#0072B2",
  "AdaPT"       = "#009E73"   
)


shapes <- c(
  "Knockoff"    = 16,  # 实心圆
  "AnnoKn-lite" = 1,   # 空心圆
  "AnnoKn"      = 8,   # 星形
  "AdaKn(RF)"   = 17,  # 实心三角
  "AdaKn(EM)"   = 15,   # 实心方形,
  "AdaPT"       = 16 
)



# 保证因子水平一致
methods_vec <- c("AdaPT","Knockoff","AdaKn(RF)", "AdaKn(EM)", "AnnoKn", "AnnoKn-lite")
plot_data$method <- factor(plot_data$method, levels = methods_vec)

# 正确的映射：shape=method（不是 shapes）
p1 <- ggplot(plot_data, aes(x = alpha, y = power, color = method, shape = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = cols, name = "Method") +
  scale_shape_manual(values = shapes, name = "Method") +
  labs(x = "Target FDR", y = "Power") +
  theme_minimal()

p2 <- ggplot(plot_data, aes(x = alpha, y = fdr, color = method, shape = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = cols, name = "Method") +
  scale_shape_manual(values = shapes, name = "Method") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  labs(x = "Target FDR", y = "FDR") +
  theme_minimal()



#library(gridExtra)
#grid.arrange(p1, p2, ncol = 2)
#combined <- grid.arrange(p1, p2, ncol = 2)


library(cowplot)

leg  <- get_legend(p1 + theme(legend.position = "right"))
p1n  <- p1 + theme(legend.position = "none")
p2n  <- p2 + theme(legend.position = "none")
prow <- plot_grid(p1n, p2n, ncol = 2, align = "hv")
final <- plot_grid(prow, leg, ncol = 2, rel_widths = c(1, 0.18))
print(final)


ggsave(
  filename = "/home/cl2667/Anno/plots/1d.png",    
  plot     = combined,         
  width    = 12,               # width
  height   = 6,                # highth
  units    = "in",             # "in"、"cm"、"mm"
  dpi      = 600               #  300 dpi
)

