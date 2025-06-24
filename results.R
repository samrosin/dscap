setting <- 1
n_t <- 3000
n_sim = 1000

wd <- "/ADD/PATH/OF/WORKING/DIRECTORY/HERE/" 
setwd(paste(wd, "setting",setting,"_nt",n_t,sep=""))
library(tidyverse)

truevalues <- read.csv(paste(wd, "truevalues_setting",
                             setting,"_nt2e+07",sep=""))
truth_rho_p <- truevalues$rho_p
truth_beta <- truevalues$beta
truth_rho_s <- truevalues$rho_s

vcov_files <- list.files()[str_detect(list.files(), "vcov")][1:n_sim]
point_files <- list.files()[!str_detect(list.files(), "vcov")][1:n_sim]

#read vcovs into list
vcov_list <- vector(mode= "list", length= length(vcov_files))
for(i in 1:length(vcov_files)){
  print(i)
  vcov_list[[i]] <- readRDS(vcov_files[i])
}

t1 = Sys.time()
# read point estimates into dataframe
point_ests <- data.frame(
  tau1_naive = rep(NA, length(point_files)),
  tau2_naive = rep(NA, length(point_files)),
  tau3_naive = rep(NA, length(point_files)),
  tau4_naive = rep(NA, length(point_files)),
  tau5_naive = rep(NA, length(point_files)),
  
  delta1_naive = rep(NA, length(point_files)),
  delta2_naive = rep(NA, length(point_files)),
  delta3_naive = rep(NA, length(point_files)),
  delta4_naive = rep(NA, length(point_files)),
  delta5_naive = rep(NA, length(point_files)),
  
  tau1 = rep(NA, length(point_files)),
  tau2 = rep(NA, length(point_files)),
  tau3 = rep(NA, length(point_files)),
  tau4 = rep(NA, length(point_files)),
  tau5 = rep(NA, length(point_files)),
  
  delta1 = rep(NA, length(point_files)),
  delta2 = rep(NA, length(point_files)),
  delta3 = rep(NA, length(point_files)),
  delta4 = rep(NA, length(point_files)),
  delta5 = rep(NA, length(point_files)),
  rho_p = rep(NA, length(point_files)),
  beta = rep(NA, length(point_files)),
  rho_s = rep(NA, length(point_files)),
  sd_rho_s = rep(NA, length(point_files)),
  rho_s_lower_percentile = rep(NA, length(point_files)),
  rho_s_upper_percentile = rep(NA, length(point_files)),
  rho_p_naive = rep(NA, length(point_files)),
  rho_s_naive = rep(NA, length(point_files)),
  beta_naive = rep(NA, length(point_files))
)

for(i in 1:length(point_files)){
  print(i)
  point_ests[i,] <- read.csv(point_files[i])
}

# bias
mean(point_ests$rho_p_naive) - truth_rho_p

mean(point_ests$rho_p) - truth_rho_p

mean(point_ests$rho_s_naive) - truth_rho_s
median(point_ests$rho_s_naive) - truth_rho_s

mean(point_ests$rho_s) - truth_rho_s
median(point_ests$rho_s) - truth_rho_s

mean(point_ests$beta_naive) - truth_beta

mean(point_ests$beta) - truth_beta

# ESEs vs. ASEs: read in ASEs from vcov files
ses_beta = rep(NA, length(vcov_list))
ses_rho_p = rep(NA, length(vcov_list))

for(i in 1:length(vcov_list)){
  ses_beta[i] = sqrt(vcov_list[[i]][53, 53])
  ses_rho_p[i] = sqrt(vcov_list[[i]][54, 54])
}

# ase vs ese: DSCAPs
sd(point_ests$beta); mean(ses_beta)
sd(point_ests$rho_p); mean(ses_rho_p)

# confidence interval coverage
point_ests$se_beta = ses_beta
point_ests$se_rho_p = ses_rho_p

point_ests = point_ests %>%
  mutate(beta_lower = beta - qnorm(.975)*se_beta,
         beta_upper = beta + qnorm(.975)*se_beta,
         rho_p_lower = rho_p - qnorm(.975)*se_rho_p,
         rho_p_upper = rho_p + qnorm(.975)*se_rho_p) %>%
  mutate(covers_beta = (beta_lower <= truth_beta) & (beta_upper >= truth_beta),
         covers_rho_p = (rho_p_lower <= truth_rho_p) & (rho_p_upper >= truth_rho_p))
mean(point_ests$covers_beta) 
mean(point_ests$covers_rho_p)




