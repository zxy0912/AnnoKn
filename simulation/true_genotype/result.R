
N.effect = 10000 # sample size
n1 = 1000 # sample size used for clustering
# p = 200 # number of variable
h2e = 0.1 # heritability 0.02, 0.05, 0.1, 0.2
alphae = 0.2 # proportion of causal SNPs. 0.05, 0.1, 0.2
p1 = 100 # number of risk SNPs
# rho = 0.5 # maximum correlation between SNPs
M = 1 # number of knockoff copies
chrid = '1'

power1 <- power2 <- power3 <- power4 <- power5 <- power6 <- numeric()
fdr1 <- fdr2 <- fdr3 <- fdr4 <- fdr5 <- fdr6 <- numeric()

for (iter in 1:100) {
  
  path = paste0("/gpfs/gibbs/pi/zhao/xz527/knockoff_anno/ghostknockoff/simulation/UKBB_simu/result/result_iter_", iter, "_n_", N.effect, "_h2e_", h2e, "_alphae_", alphae, "_p1_", p1, "_M_", M, "_chr_", chrid, ".RData")
  load(path)
  
  power1 <- rbind(power1, result$power1)
  power2 <- rbind(power2, result$power2)
  power3 <- rbind(power3, result$power3)
  power4 <- rbind(power4, result$power4)
  power5 <- rbind(power5, result$power5)
  power6 <- rbind(power6, result$power6)
  
  fdr1 <- rbind(fdr1, result$fdr1)
  fdr2 <- rbind(fdr2, result$fdr2)
  fdr3 <- rbind(fdr3, result$fdr3)
  fdr4 <- rbind(fdr4, result$fdr4)
  fdr5 <- rbind(fdr5, result$fdr5)
  fdr6 <- rbind(fdr6, result$fdr6)
}

colMeans(power1, na.rm = TRUE)
colMeans(power2, na.rm = TRUE)
colMeans(power3, na.rm = TRUE)

colMeans(fdr1, na.rm = TRUE)
colMeans(fdr2, na.rm = TRUE)
colMeans(fdr3, na.rm = TRUE)

colMeans(power4, na.rm = TRUE)
colMeans(power5, na.rm = TRUE)
colMeans(power6, na.rm = TRUE)


colMeans(fdr4, na.rm = TRUE)
colMeans(fdr5, na.rm = TRUE)
colMeans(fdr6, na.rm = TRUE)








library(ggplot2)

methods <- c("knockoff", "knockoff-simple", "knockoff-anno", 
             "GhostKnockoff", "GK-simple", "GK-anno")

fdr_seq <- seq(0.4, 0.1, by = -0.05)

power_means <- list(
  colMeans(power1[,1:7], na.rm = TRUE),
  colMeans(power2[,1:7], na.rm = TRUE),
  colMeans(power3[,1:7], na.rm = TRUE),
  colMeans(power4[,1:7], na.rm = TRUE),
  colMeans(power5[,1:7], na.rm = TRUE),
  colMeans(power6[,1:7], na.rm = TRUE)
)

fdr_means <- list(
  colMeans(fdr1[,1:7], na.rm = TRUE),
  colMeans(fdr2[,1:7], na.rm = TRUE),
  colMeans(fdr3[,1:7], na.rm = TRUE),
  colMeans(fdr4[,1:7], na.rm = TRUE),
  colMeans(fdr5[,1:7], na.rm = TRUE),
  colMeans(fdr6[,1:7], na.rm = TRUE)
)

df_power <- do.call(rbind, lapply(1:6, function(i) {
  data.frame(Method = methods[i],
             q_value = fdr_seq,
             Value = power_means[[i]],
             Metric = "Power")
}))

df_fdr <- do.call(rbind, lapply(1:6, function(i) {
  data.frame(Method = methods[i],
             q_value = fdr_seq,
             Value = fdr_means[[i]],
             Metric = "FDR")
}))

df <- rbind(df_power, df_fdr)


df$Metric <- factor(df$Metric, levels = c("Power", "FDR"))

# ggplot(df, aes(x = q_value, y = Value, color = Method)) +
#   geom_line(size = 1) +
#   geom_point(size = 2) +
#   facet_wrap(~Metric, ncol = 2, scales = "free_y") +
#   geom_line(data = data.frame(q_value = fdr_seq, Value = fdr_seq, Metric = factor("FDR", levels = c("Power","FDR"))),
#             aes(x = q_value, y = Value),
#             inherit.aes = FALSE, color = "black", linetype = "dashed") +
#   theme_minimal(base_size = 14) +
#   labs(x = "q-value threshold", y = "Mean Value", 
#        title = paste0("n_", N.effect, "_h2e_", h2e, "_alphae_", alphae, "_p1_", p1)) +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 30, hjust = 1))



ggplot(df, aes(x = q_value, y = Value, 
               color = Method, shape = Method, linetype = Method)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~Metric, ncol = 2, scales = "free_y") +
  geom_line(
    data = data.frame(
      q_value = fdr_seq, 
      Value = fdr_seq, 
      Metric = factor("FDR", levels = c("Power","FDR"))
    ),
    aes(x = q_value, y = Value),
    inherit.aes = FALSE, 
    color = "black", linetype = "dashed"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    x = "q-value threshold", 
    y = "Mean Value", 
    title = paste0("n_", N.effect, "_h2_", h2e, "_alphae_", alphae, "_p1_", p1)
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

