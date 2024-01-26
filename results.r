setting <- 3
n_t <- 3000
epsilon <- .0000001
setwd(paste("/nas/longleaf/home/srosin/paper2/setting",setting,"_nt",n_t,sep=""))
n_sim = length(list.files())

library(tidyverse)

truevalues <- read.csv(paste("/nas/longleaf/home/srosin/paper2/truevalues/truevalues_setting",
                  setting,"_nt2e+07",sep=""))

truth_rho_p <- truevalues$rho_p
truth_rho_s <- truevalues$rho_s
truth_beta <- truevalues$beta

results <- data.frame(rho_p = rep(NA, n_sim),
                      #sd_rho_p = rep(NA, n_sim),
                      #rho_p_lower = rep(NA, n_sim),
                      #rho_p_upper = rep(NA, n_sim),
                      beta = rep(NA, n_sim),
                     # sd_beta = rep(NA, n_sim),
                     # beta_lower = rep(NA, n_sim),
                     # beta_upper = rep(NA, n_sim),
                      rho_s = rep(NA, n_sim),
                      sd_rho_s = rep(NA, n_sim),
                      rho_s_lower_percentile = rep(NA, n_sim),
                      rho_s_upper_percentile = rep(NA, n_sim),
                      rho_s_lower_bca = rep(NA, n_sim),
                      rho_s_upper_bca = rep(NA, n_sim),
                     rho_s_lower_bca_coxed = rep(NA, n_sim),
                     rho_s_upper_bca_coxed = rep(NA, n_sim),
                      rho_p_naive = rep(NA, n_sim),
                      rho_s_naive = rep(NA, n_sim),
                      beta_naive = rep(NA, n_sim))

i <- 1
for(i in 1:length(list.files())){
  results[i,] <- read.csv(list.files()[i])
  i <- i + 1
}

##### ONLY KEEP 100 
results <- results[1:100,]

results <- results %>% 
  mutate(covers_beta = (results$beta_lower <= truth_beta + epsilon) & 
           (results$beta_upper >= truth_beta - epsilon),
         covers_rho_p = (results$rho_p_lower <= truth_rho_p + epsilon) & 
           (results$rho_p_upper >= truth_rho_p + epsilon),
    covers_rho_s = (results$rho_s_lower_percentile <= truth_rho_s + epsilon) & 
           (results$rho_s_upper_percentile >= truth_rho_s - epsilon),
         covers_rho_s_bca = (results$rho_s_lower_bca <= truth_rho_s + epsilon) & 
           (results$rho_s_upper_bca >= truth_rho_s - epsilon))

mean(results$beta) - truth_beta
mean(results$rho_p) - truth_rho_p
mean(results$rho_s) - truth_rho_s

median(results$beta) - truth_beta
median(results$rho_p) - truth_rho_p
median(results$rho_s) - truth_rho_s

# histograms -- are the distributions approximately normal? can use shapiro wilks test
hist(results$beta, breaks = 20)
hist(results$rho_p, breaks = 20)
hist(results$rho_s)

sd(results$beta); mean(results$sd_beta)
sd(results$rho_p); mean(results$sd_rho_p)
sd(results$rho_s); mean(results$sd_rho_s)

mean(results$covers_beta)
mean(results$covers_rho_p)
mean(results$covers_rho_s)
mean(results$covers_rho_s_bca, na.rm = T)

mean(results$beta_naive) - truth_beta
mean(results$rho_p_naive) - truth_rho_p
mean(results$rho_s_naive) - truth_rho_s

## plot a sample of sims with error bars

results_ordered <- results[sample(nrow(results), 200),]
results_ordered <- results_ordered[order(results_ordered$rho_s),] %>% 
  mutate(sample = 1:200)

rho_s_plot <- results_ordered %>% 
  ggplot(aes(x = sample, y = rho_s)) + 
  geom_point(aes(color = covers_rho_s_bca)) + 
  geom_errorbar(aes(ymin = rho_s_lower_bca, ymax = rho_s_upper_bca, 
                    color = covers_rho_s), alpha = .4) + 
  coord_flip() + 
  geom_hline(yintercept = .01, linetype = "dashed", color = "black") + 
  #labs(title = "DGP 4: 95% Confidence Intervals") + theme(plot.title = element_text(hjust = 0.5)) + ylim(0, .4) +  
  scale_color_manual(values = c("black", "gray"), labels = c("No", "Yes")) + 
  labs(colour="CI covers?") +
  xlab("Ordered simulation index") + 
  ylab(expression(hat(rho)[s])) + 
  theme_classic() + 
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = .5),
        axis.title = element_text(size = 20),
        legend.position = c(.9, 0.8),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 16))
rho_s_plot

