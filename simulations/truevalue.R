wd <- "/ADD/PATH/OF/WORKING/DIRECTORY/HERE/" 
setwd(wd)
library(tidyverse)
library(MASS)

set.seed(123)
num_trials <- 5

#collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
setting <- as.numeric(args[1])
n_t <- as.numeric(args[2])

source(paste(wd, "parameter_values.R",sep=""))

### outcome models
Ymod = Y ~ X1 + X2 
Smod = S ~ X1 + X2 

inv_logit <- function(x){exp(x) / (1 + exp(x))}

# generate data, only necessary to generate for trial 1 (which is a random sample from the target population)-----------------------------------------------------------

# Trial variable - for fixed number of trials 
trial <- rep(1, n_t)

# random treatment assignment -- could vary by trial
A <- rbinom(n=n_t, size = 1, prob = 0.5)
df <- data.frame(trial, A) %>%
  rowwise() %>%
  mutate(A = ifelse(A==0, 0, trial))

# target population indicator
df$R <- ifelse(trial == 1, 0, 1)

# covariate generation
df$X1 <- NA; df$X2 <- NA
df$X1 <- rbinom(n = nrow(df), size = 1, prob = 0.3 )

df$X2 <- rnorm(n_t, mean = 55, sd = 22.3)

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

ey0 <- mean(df$Y_0)
ey1 <- mean(df$Y_1)
ey2 <- mean(df$Y_2)
ey3 <- mean(df$Y_3)
ey4 <- mean(df$Y_4)
ey5 <- mean(df$Y_5)

es0 <- mean(df$S_0)
es1 <- mean(df$S_1)
es2 <- mean(df$S_2)
es3 <- mean(df$S_3)
es4 <- mean(df$S_4)
es5 <- mean(df$S_5)

tau <- c(1 - ey1 / ey0, 1 - ey2 / ey0, 1 - ey3 / ey0, 1 - ey4 / ey0, 1 - ey5 / ey0)
delta <- c(es1 - es0, es2 - es0, es3 - es0, es4 - es0, es5 - es0)

rho_p <- cor(tau, delta, method = "pearson")
rho_s <- cor(tau, delta, method = "spearman")
beta <- coef(lm(tau ~ delta))[2]

truevalues <- data.frame(n_t, rho_p, rho_s, beta,
                         tau[1], tau[2], tau[3], tau[4], tau[5],
                         delta[1], delta[2], delta[3], delta[4], delta[5],
                         ey0, ey1, ey2, ey3, ey4, ey5, 
                         es0, es1, es2, es3, es4, es5)

output_filename <- paste(wd,"truevalues/truevalues_setting",setting,"_nt",n_t,sep="")
write.csv(truevalues, output_filename)
