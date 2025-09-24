library(ggplot2)

iteration = 100

alphalist <- seq(0.4, 0.1, by = -0.05)
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
power9 = matrix(0, nrow = iteration, ncol = len)
power10 = matrix(0, nrow = iteration, ncol = len)
power11 = matrix(0, nrow = iteration, ncol = len)

power_em = matrix(0, nrow = iteration, ncol = 1)
power_em_gk = matrix(0, nrow = iteration, ncol = 1)

fdr1 = matrix(0, nrow = iteration, ncol = len)
fdr2 = matrix(0, nrow = iteration, ncol = len)
fdr3 = matrix(0, nrow = iteration, ncol = len)
fdr4 = matrix(0, nrow = iteration, ncol = len)
fdr5 = matrix(0, nrow = iteration, ncol = len)
fdr6 = matrix(0, nrow = iteration, ncol = len)
fdr7 = matrix(0, nrow = iteration, ncol = len)
fdr8 = matrix(0, nrow = iteration, ncol = len)
fdr9 = matrix(0, nrow = iteration, ncol = len)
fdr10 = matrix(0, nrow = iteration, ncol = len)
fdr11 = matrix(0, nrow = iteration, ncol = len)
fdr_em = matrix(0, nrow = iteration, ncol = 1)
fdr_em_gk =matrix(0, nrow = iteration, ncol = 1)


s = 0
removelist = c()

binary = ''
dimen = 'double'
heri = 0.05
N.effect = 5000 #
p = 300

for (i in 1:iteration) {
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/ghostknockoff/simulation/AnnoGK_simu/result/heri_",heri,"_n_",N.effect,"_p_",p, "_", i,".RData")
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
  power10[i,] <- result$power10
  power11[i,] <- result$power11
  
  fdr1[i,]   <- result$fdr1
  fdr2[i,]   <- result$fdr2
  fdr3[i,]   <- result$fdr3
  fdr4[i,]   <- result$fdr4
  fdr5[i,]   <- result$fdr5
  fdr6[i,]   <- result$fdr6
  fdr7[i,]   <- result$fdr7
  fdr8[i,]   <- result$fdr8
  fdr9[i,]   <- result$fdr9
  fdr10[i,]   <- result$fdr10
  fdr11[i,]   <- result$fdr11
  
}

# if(length(removelist) > 0){
#   power1 = power1[-removelist,]
#   power2 = power2[-removelist,]
#   power3 = power3[-removelist,]
#   power4 = power4[-removelist,]
#   power5 = power5[-removelist,]
#   power6 = power6[-removelist,]
#   power7 = power7[-removelist,]
#   power8 = power8[-removelist,]
#   fdr1 = fdr1[-removelist,]
#   fdr2 = fdr2[-removelist,]
#   fdr3 = fdr3[-removelist,]
#   fdr4 = fdr4[-removelist,]
#   fdr5 = fdr5[-removelist,]
#   fdr6 = fdr6[-removelist,]
#   fdr7 = fdr7[-removelist,]
#   fdr8 = fdr8[-removelist,]
# }


power <- c(
  colMeans(power1, na.rm = TRUE),
  # colMeans(power2, na.rm = TRUE),
  colMeans(power3, na.rm = TRUE),
  colMeans(power4, na.rm = TRUE),  
  # colMeans(power5, na.rm = TRUE),
  colMeans(power6, na.rm = TRUE),  
  colMeans(power7, na.rm = TRUE),
  colMeans(power8, na.rm = TRUE),
  # colMeans(power9, na.rm = TRUE), 
  colMeans(power10, na.rm = TRUE),
  colMeans(power11, na.rm = TRUE) 
)


fdr <- c(
  colMeans(fdr1, na.rm = TRUE),
  # colMeans(fdr2, na.rm = TRUE),
  colMeans(fdr3, na.rm = TRUE),
  colMeans(fdr4, na.rm = TRUE),
  # colMeans(fdr5, na.rm = TRUE),
  colMeans(fdr6, na.rm = TRUE),
  colMeans(fdr7, na.rm = TRUE),
  colMeans(fdr8, na.rm = TRUE),
  # colMeans(fdr9, na.rm = TRUE),
  colMeans(fdr10, na.rm = TRUE),
  colMeans(fdr11, na.rm = TRUE)
)


name_method <- c('Knockoff', 
                 # 'AnnoKn-simple', 
                 'AnnoKn', 
                 'GhostKnockoff', 
                 # 'AnnoGK-simple',
                 'AnnoGK',
                 'AnnoGK-dss', 
                 'GhostKnockoff M=5', 
                 # 'AnnoGK-simple M=5',
                 'AnnoGK M=5', 
                 'AnnoGK-dss M=5')

method <- rep(name_method, each = len)
plot_data <- data.frame(alpha = rep(alphalist, 8),
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

