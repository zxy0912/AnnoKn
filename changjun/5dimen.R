library(ggplot2)

iteration = 100

alphalist <- seq(0.3,0.05,-0.05)
len = length(alphalist)

power1 = matrix(0, nrow = iteration, ncol = len)
power2 = matrix(0, nrow = iteration, ncol = len)
power3 = matrix(0, nrow = iteration, ncol = len)
power4 = matrix(0, nrow = iteration, ncol = len)
power5 = matrix(0, nrow = iteration, ncol = len)
power6 = matrix(0, nrow = iteration, ncol = len)
power7 = matrix(0, nrow = iteration, ncol = len)
power8 = matrix(0, nrow = iteration, ncol = len)
power9 = matrix(0, nrow = iteration, ncol = len)
lambda_inf = matrix(0, nrow = iteration, ncol = 2)
lambda_noninf = matrix(0, nrow = iteration, ncol = 8)

fdr1 = matrix(0, nrow = iteration, ncol = len)
fdr2 = matrix(0, nrow = iteration, ncol = len)
fdr3 = matrix(0, nrow = iteration, ncol = len)
fdr4 = matrix(0, nrow = iteration, ncol = len)
fdr5 = matrix(0, nrow = iteration, ncol = len)
fdr6 = matrix(0, nrow = iteration, ncol = len)
fdr7 = matrix(0, nrow = iteration, ncol = len)
fdr8 = matrix(0, nrow = iteration, ncol = len)
fdr9 = matrix(0, nrow = iteration, ncol = len)
fdr_em = matrix(0, nrow = iteration, ncol = 1)

s = 0
removelist = c()

binary = 'binary'
dimen = 'double'
amp = 25
type = 'in'

for (i in 1:iteration) {
  path = paste0("/gpfs/gibbs/project/zhao/cl2667/knockoff_side/iter_simu/double/Sen/10dimS/result_", binary, type, "_amp_",amp, "_", i,".RData")
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
  power9[i,] <- result$power9
  lambda_inf[i,] <- result$lambda_inf
  lambda_noninf[i,] <- result$lambda_noninf
  
  
  fdr1[i,]   <- result$fdr1
  fdr2[i,]   <- result$fdr2
  fdr3[i,]   <- result$fdr3
  fdr4[i,]   <- result$fdr4
  fdr5[i,]   <- result$fdr5
  fdr6[i,]   <- result$fdr6
  fdr7[i,]   <- result$fdr7
  fdr8[i,]   <- result$fdr8
  fdr9[i,]   <- result$fdr9
  
}


if(length(removelist) > 0){
  power1 = power1[-removelist,]
  power2 = power2[-removelist,]
  power3 = power3[-removelist,]
  power4 = power4[-removelist,]
  power5 = power5[-removelist,]
  power6 = power6[-removelist,]
  power7 = power7[-removelist,]
  power8 = power8[-removelist,]
  power9 = power9[-removelist,]
  lambda_inf = lambda_inf[-removelist,]
  lambda_noninf = lambda_noninf[-removelist,]
  
  fdr1 = fdr1[-removelist,]
  fdr2 = fdr2[-removelist,]
  fdr3 = fdr3[-removelist,]
  fdr4 = fdr4[-removelist,]
  fdr5 = fdr5[-removelist,]
  fdr6 = fdr6[-removelist,]
  fdr7 = fdr7[-removelist,]
  fdr8 = fdr8[-removelist,]
  fdr9 = fdr9[-removelist,]
}


power <- c(
  colMeans(power1, na.rm = TRUE), # "Knockoff"
  #colMeans(power2, na.rm = TRUE),
  colMeans(power3, na.rm = TRUE),
  #colMeans(power4, na.rm = TRUE),
  colMeans(power5, na.rm = TRUE),
  #colMeans(power6, na.rm = TRUE),
  colMeans(power7, na.rm = TRUE),
  colMeans(power9, na.rm = TRUE)
)


fdr <- c(
  colMeans(fdr1, na.rm = TRUE),
  #colMeans(fdr2, na.rm = TRUE),
  colMeans(fdr3, na.rm = TRUE),
  #colMeans(fdr4, na.rm = TRUE),
  colMeans(fdr5, na.rm = TRUE),
  #colMeans(fdr6, na.rm = TRUE),
  colMeans(fdr7, na.rm = TRUE),
  colMeans(fdr9, na.rm = TRUE)
)

method <- rep(c("Knockoff",  "AnnoKn with mixed annotation", "AnnoKn with informative annotation" ,"AdaKn(EM) with mixed annotation","AdaKn(EM) with informative annotation"), each = len)
plot_data <- data.frame(alpha = rep(alphalist, 5),
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

library(gridExtra)
grid.arrange(p1, p2, ncol = 2)


combined <- grid.arrange(p1, p2, ncol = 2)

cols <- c(
  "Knockoff"                               = "#E69F00",
  "AdaKn(EM) with informative annotation"  = "#0072B2", # 深蓝
  "AdaKn(EM) with mixed annotation"        = "#56B4E9", # 浅蓝
  "AdaKn(RF) with informative annotation"  = "#D55E00", # 深橙
  "AdaKn(RF) with mixed annotation"        = "#F4A582", # 浅橙
  "AnnoKn with informative annotation"     = "#FF0000", # 深红
  "AnnoKn with mixed annotation"           = "#FB8072"  # 浅红/鲑色
)

shapes <- c(
  "Knockoff"                               = 16, # 实心圆
  "AdaKn(EM) with informative annotation"  = 15, # 实心方
  "AdaKn(EM) with mixed annotation"        = 0,  # 空心方
  "AdaKn(RF) with informative annotation"  = 17, # 实心三角
  "AdaKn(RF) with mixed annotation"        = 2,  # 空心三角
  "AnnoKn with informative annotation"     = 8,  # 星形
  "AnnoKn with mixed annotation"           = 5   # 菱形
)

methods_vec <- c("AdaKn(EM) with mixed annotation","AdaKn(EM) with informative annotation", "AnnoKn with mixed annotation", "AnnoKn with informative annotation","Knockoff")
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

library(cowplot)



library(cowplot)

leg <- get_legend(p1 + theme(legend.position = "right"))

p1n <- p1 + theme(legend.position = "none")
p2n <- p2 + theme(legend.position = "none")

row1 <- plot_grid(p1n, ncol = 1, align = "hv")
row2 <- plot_grid(p2n, ncol = 1, align = "hv")

plots <- plot_grid(row1, row2, ncol = 1)

final <- plot_grid(plots, leg, ncol = 2, rel_widths = c(0.6, 0.4))

print(final)








library(ggplot2)

# 压平成向量
lambda_noninf_vec <- as.numeric(lambda_noninf)
lambda_inf_vec    <- as.numeric(lambda_inf)

# 组合到同一数据框
df <- rbind(
  data.frame(value = lambda_noninf_vec, group = "non-informative"),
  data.frame(value = lambda_inf_vec,    group = "informative")
)

# 绘图
ggplot(na.omit(df), aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.5,
                 breaks = seq(-5, 15, by = 0.25),
                 boundary = 0) +
  scale_fill_manual(values = c("non-informative" = "#E69F00",
                               "informative"     = "#0072B2")) +
  labs(
    x = expression(lambda[l]),  # ✅ 数学公式标签
    y = "Frequency",
    fill = "Group"
  ) +
  coord_cartesian(xlim = c(-5, 15)) +
  theme_minimal()


