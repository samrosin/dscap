iter <- Sys.getenv("SLURM_ARRAY_TASK_ID") # get the array ID to use as the random seed
set.seed(iter)
#sim=1
#set.seed(123)
library(tidyverse)
library(geex)
library(MASS)
library(lmtest)
library(coxed)
wd <- "/nas/longleaf/home/srosin/paper2/"
setwd(wd)

#collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
setting <- as.numeric(args[1])
n_t <- as.numeric(args[2]) # number of participants per trial
n_boot <- as.numeric(args[3]) # number of bootstrap samples

source(paste(wd, "estFUN_taudelta.R",sep=""))
source(paste(wd, "parameter_values.R",sep=""))

# more parameters /arguments; could pass these from slurm
num_trials <- 5
alpha <- .05 # alpha level for confidence intervals
deriv_method <- "simple" # "simple" or "richardson"


### outcome models
Ymod = Y ~ X1 + X2 
Smod = S ~ X1 + X2 

inv_logit <- function(x){exp(x) / (1 + exp(x))}

# generate data -----------------------------------------------------------

# Trial variable - for fixed number of trials 5
trial <- c(rep(1, n_t), rep(2, n_t), rep(3, n_t), 
           rep(4, n_t), rep(5, n_t))

# random treatment assignment -- could vary by trial
A <- rbinom(n=n_t*num_trials, size = 1, prob = 0.5)
df <- data.frame(trial, A) %>%
  rowwise() %>%
  mutate(A = ifelse(A==0, 0, trial))

# target population indicator
df$R <- ifelse(trial == 1, 0, 1)

# covariate generation
df$X1 <- NA; df$X2 <- NA
df$X1 <- rbinom(n = nrow(df), size = 1,
                prob = 0.3 - 0.1*(df$trial==2) - 0.1*(df$trial==3)
                + 0.1*(df$trial==4) + 0.1*(df$trial==5))

df[df$trial==1,"X2"] <- rnorm(n_t, mean = 55, sd = 22.3)
df[(df$trial==2) | (df$trial==4),"X2"] <- rnorm(n_t*2, mean = 45, sd = 18)
df[(df$trial==3) | (df$trial==5),"X2"] <- rnorm(n_t*2, mean = 65, sd = 24)

# truncate age (X2) into [18, 85]
df[df$X2<18,]$X2 <- 18
df[df$X2>85,]$X2 <- 85

# create covariate matrix for potential outcome generation
intercept <- matrix(rep(1, nrow(df)),ncol=1)
X <- as.matrix(df[,c("X1","X2")])
X <- cbind(intercept, X)

# generate potential outcomes
for(a in 1:(num_trials+1)){
  Y_a <- paste("Y_", a-1, sep="")
  S_a <- paste("S_", a-1, sep="")
  theta_a <- matrix(theta[[a]], ncol=1)
  gamma_a <- matrix(gamma[[a]], ncol = 1)
  df[[Y_a]] <- rbinom(nrow(df), 1, prob = inv_logit(X %*% theta_a))
  if(a == 1){
    df[[S_a]] <- rnorm(nrow(df), mean = (X %*% gamma_a), sd = 0.28)
  } else {
    df[[S_a]] <- rnorm(nrow(df), mean = (X %*% gamma_a), sd = 0.64)
  }
}

# generate outcomes by consistency
df$Y <- NA
df$S <- NA
for(a in 1:(num_trials+1)){
  Y_a <- paste("Y_", a-1, sep="")
  df[df$A==(a-1),]$Y <- unlist(df[df$A==(a-1),Y_a])
  S_a <- paste("S_", a-1, sep="")
  df[df$A==(a-1),]$S <- unlist(df[df$A==(a-1),S_a])
}

# keep necessary variables
df <- df %>% dplyr::select(trial, A, R, X1, X2, Y, S)
df_r1 <- df %>% filter(R == 1)

# outcome models E[O | X=x, A=a]
model_y0 <- glm(Ymod, data=df, family = binomial, subset = (A==0))
model_y2 <- glm(Ymod, data=df, family = binomial, subset = (A==2))
model_y3 <- glm(Ymod, data=df, family = binomial, subset = (A==3))
model_y4 <- glm(Ymod, data=df, family = binomial, subset = (A==4))
model_y5 <- glm(Ymod, data=df, family = binomial, subset = (A==5))

model_s0 <- glm(Smod, data=df, family = gaussian, subset = (A==0))
model_s2 <- glm(Smod, data=df, family = gaussian, subset = (A==2))
model_s3 <- glm(Smod, data=df, family = gaussian, subset = (A==3))
model_s4 <- glm(Smod, data=df, family = gaussian, subset = (A==4))
model_s5 <- glm(Smod, data=df, family = gaussian, subset = (A==5))

df$py0 <- predict(model_y0,newdata=df, type="response")
df$py2 <- predict(model_y2,newdata=df, type="response")
df$py3 <- predict(model_y3,newdata=df, type="response")
df$py4 <- predict(model_y4,newdata=df, type="response")
df$py5 <- predict(model_y5,newdata=df, type="response")

df$ps0 <- predict(model_s0,newdata=df, type="response")
df$ps2 <- predict(model_s2,newdata=df, type="response")
df$ps3 <- predict(model_s3,newdata=df, type="response")
df$ps4 <- predict(model_s4,newdata=df, type="response")
df$ps5 <- predict(model_s5,newdata=df, type="response")

# naive estimates
ey0_naive <- mean(df[(df$A == 0),]$Y)
ey1_naive <- mean(df[(df$A == 1),]$Y)
ey2_naive <- mean(df[(df$A == 2),]$Y)
ey3_naive <- mean(df[(df$A == 3),]$Y)
ey4_naive <- mean(df[(df$A == 4),]$Y)
ey5_naive <- mean(df[(df$A == 5),]$Y)

es0_naive <- mean(df[(df$A == 0),]$S)
es1_naive <- mean(df[(df$A == 1),]$S)
es2_naive <- mean(df[(df$A == 2),]$S)
es3_naive <- mean(df[(df$A == 3),]$S)
es4_naive <- mean(df[(df$A == 4),]$S)
es5_naive <- mean(df[(df$A == 5),]$S)

# estimate E[Y^a | R = 0] and E[S^a | R = 0]
ey1 <- sum(df[df$A==1,]$Y) / (nrow(df[df$A==1,]))
es1 <- sum(df[df$A==1,]$S) / (nrow(df[df$A==1,]))

n_r <- sum(1 - df$R)
ey0 <- (1 / n_r) * sum((1-df$R)*df$py0 )
ey2 <- (1 / n_r) * sum((1-df$R)*df$py2 )
ey3 <- (1 / n_r) * sum((1-df$R)*df$py3 )
ey4 <- (1 / n_r) * sum((1-df$R)*df$py4 )
ey5 <- (1 / n_r) * sum((1-df$R)*df$py5 )

es0 <- (1 / n_r) * sum((1-df$R)*df$ps0 )
es2 <- (1 / n_r) * sum((1-df$R)*df$ps2 )
es3 <- (1 / n_r) * sum((1-df$R)*df$ps3 )
es4 <- (1 / n_r) * sum((1-df$R)*df$ps4 )
es5 <- (1 / n_r) * sum((1-df$R)*df$ps5 )

tau_naive <- c(1 - ey1_naive / ey0_naive, 1 - ey2_naive / ey0_naive, 1 - ey3_naive / ey0_naive, 
               1 - ey4_naive / ey0_naive, 1 - ey5_naive / ey0_naive)
delta_naive <- c(es1_naive - es0_naive, es2_naive - es0_naive, es3_naive - es0_naive, 
                 es4_naive - es0_naive, es5_naive - es0_naive)
rho_p_naive <- cor(tau_naive, delta_naive, method = "pearson")
rho_s_naive <- cor(tau_naive, delta_naive, method = "spearman")
beta_naive <- coef(lm(tau_naive ~ delta_naive))[2]

tau <- c(1 - ey1 / ey0, 1 - ey2 / ey0, 1 - ey3 / ey0, 1 - ey4 / ey0, 1 - ey5 / ey0)
delta <- c(es1 - es0, es2 - es0, es3 - es0, es4 - es0, es5 - es0)
rho_p <- cor(tau, delta, method = "pearson")
rho_s <- cor(tau, delta, method = "spearman")
beta <- coef(lm(tau ~ delta))[2]

# store models in a list
models <- list(model_y0 = model_y0,
               model_y2 = model_y2, model_y3 = model_y3,
               model_y4 = model_y4, model_y5 = model_y5,
               model_s0 = model_s0,
               model_s2 = model_s2, model_s3 = model_s3,
               model_s4 = model_s4, model_s5 = model_s5)

# coefficients from regression models
coef_y0 <- coef(models$model_y0)
coef_y2 <- coef(models$model_y2)
coef_y3 <- coef(models$model_y3)
coef_y4 <- coef(models$model_y4)
coef_y5 <- coef(models$model_y5)

coef_s0 <- coef(models$model_s0)
coef_s2 <- coef(models$model_s2)
coef_s3 <- coef(models$model_s3)
coef_s4 <- coef(models$model_s4)
coef_s5 <- coef(models$model_s5)

theta <- c(coef_y0, coef_y2, coef_y3, coef_y4, coef_y5,
           coef_s0, coef_s2, coef_s3, coef_s4, coef_s5,
           ey0, ey1, ey2, ey3, ey4, ey5,
           es0, es1, es2, es3, es4, es5,
           tau[1], tau[2], tau[3], tau[4], tau[5],
           delta[1], delta[2], delta[3], delta[4], delta[5],
           rho_p, beta)

m_est <-m_estimate(
  estFUN = estFUN_taudelta,
  data  = df,
  inner_args = list(models = models),
  roots = theta,
  compute_roots = F,
  deriv_control = setup_deriv_control(method = deriv_method),
  compute_vcov = T
)

vcov_m_est <- vcov(m_est)

tau1 <- tau[1]; tau2 <- tau[2]; tau3 <- tau[3]; tau4 <- tau[4]; tau5 <- tau[5]
delta1 <- delta[1]; delta2 <- delta[2]; delta3 <- delta[3]; delta4 <- delta[4]; delta5 <- delta[5]

tau1_naive <- tau_naive[1]; tau2_naive <- tau_naive[2]; tau3_naive <- tau_naive[3]
tau4_naive <- tau_naive[4]; tau5_naive <- tau_naive[5]

delta1_naive <- delta_naive[1]; delta2_naive <- delta_naive[2]; delta3_naive <- delta_naive[3]
delta4_naive <- delta_naive[4]; delta5_naive <- delta_naive[5]

results <- data.frame(tau1_naive,tau2_naive,tau3_naive,tau4_naive,tau5_naive,
                      delta1_naive,delta2_naive,delta3_naive,delta4_naive,delta5_naive,
                      tau1,tau2,tau3,tau4,tau5,
                      delta1,delta2,delta3,delta4,delta5,
                      rho_p,
                      #sd_rho_p,
                      #rho_p_lower = rho_p - qnorm(1 - alpha/2)*sd_rho_p,
                      #rho_p_upper = rho_p + qnorm(1 - alpha/2)*sd_rho_p,
                      beta,
                     # sd_beta,
                     # beta_lower = beta - qnorm(1-alpha/2)*sd_beta,
                     # beta_upper = beta + qnorm(1-alpha/2)*sd_beta,
                      rho_s, 
                      #sd_rho_s,
                      #rho_s_lower_percentile,
                      #rho_s_upper_percentile,
                      # rho_s_lower_bca,
                      # rho_s_upper_bca,
                      rho_p_naive,
                      rho_s_naive,
                      beta_naive
)

output_dir <- paste(wd,"setting",setting,"_nt",n_t,sep="")
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

vcov_filename <- paste(output_dir,"/iter_",iter,"_vcov.RDS",sep="")
saveRDS(vcov_m_est, vcov_filename)

output_filename <- paste(output_dir,"/iter_",iter,".csv",sep="")
write_csv(results, output_filename)
