iter <- Sys.getenv("SLURM_ARRAY_TASK_ID") # get the array ID to use as the random seed
set.seed(iter)

wd <- "/ADD/PATH/OF/WORKING/DIRECTORY/HERE/" 
setwd(wd)
library(tidyverse)
library(geex)
library(MASS)
library(lmtest)

#collect arguments passed in from SLURM job script
args <- commandArgs(trailingOnly=TRUE)
setting <- as.numeric(args[1])
n_t <- as.numeric(args[2]) # number of participants per trial
n_perm <- as.numeric(args[3]) # number of permutations

if(setting %in% c(1, 2, 3, 4)){
  source(paste(wd, "parameter_values.R", sep=""))
} else if(setting == 5){
  source(paste(wd, "parameter_values_nonequivalent_placebo.R", sep=""))
}

# more parameters /arguments; could pass these from slurm
num_trials <- 5

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
if(setting %in% 1:4){
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
} else if(setting == 5){
  for(t in 1:num_trials){
    Y_0t <- paste("Y_0", t, sep="")
    theta_0t <- matrix(unlist(theta[[t]][1]), ncol = 1)
    df[[Y_0t]] <- rbinom(nrow(df), 1, prob = inv_logit(X %*% theta_0t))
    
    Y_at <- paste("Y_", t, t, sep="")
    theta_at <- matrix(unlist(theta[[t]][2]), ncol = 1)
    df[[Y_at]] <- rbinom(nrow(df), 1, prob = inv_logit(X %*% theta_at))
    
    S_0t <- paste("S_0", t, sep="")
    gamma_0t <- matrix(unlist(gamma[[t]][1]), ncol = 1)
    df[[S_0t]] <- rnorm(nrow(df), mean = (X %*% gamma_0t), sd = .64)
    
    S_at <- paste("S_", t, t, sep="")
    gamma_at <- matrix(unlist(gamma[[t]][2]), ncol = 1)
    df[[S_at]] <- rnorm(nrow(df), mean = (X %*% gamma_at), sd = .64)
  }
}

# generate outcomes by consistency
df$Y <- NA
df$S <- NA

if(setting %in% 1:4){
  for(a in 1:(num_trials+1)){
    Y_a <- paste("Y_", a-1, sep="")
    df[df$A==(a-1),]$Y <- unlist(df[df$A==(a-1),Y_a])
    S_a <- paste("S_", a-1, sep="")
    df[df$A==(a-1),]$S <- unlist(df[df$A==(a-1),S_a])
  }
} else if(setting == 5){
  for(a in unique(df$A)){
    for(t in unique(df$trial)){
      Y_at <- paste("Y_",a,t,sep="")
      S_at <- paste("S_",a,t,sep="")
      if(nrow(df[(df$A==a) & (df$trial==t),]) > 0){ # for all combinations {A=a,T=t} under consideration
        df[(df$A==a) & (df$trial==t),]$Y <- unlist(df[(df$A==a) & (df$trial==t),Y_at]) 
        df[(df$A==a) & (df$trial==t),]$S <- unlist(df[(df$A==a) & (df$trial==t),S_at]) 
      }
    }
  }
}


# keep necessary variables
df <- df %>% dplyr::select(trial, A, R, X1, X2, Y, S)

# # likelihood ratio tests
df_placebo <- df %>% filter(A == 0)

ymod_no_trial <- glm(Y ~ X1 + X2 , binomial, df_placebo)
ymod_trial <- glm(Y ~ X1 + X2 + trial + X1:trial + X2:trial, binomial, df_placebo)
lrt_y <- lmtest::lrtest(ymod_no_trial, ymod_trial)[2,4] #LRT test statistic
lrt_y_pval <- lmtest::lrtest(ymod_no_trial, ymod_trial)[2,5] #extract LRT p-value

smod_no_trial <- glm(S ~ X1 + X2 , gaussian, df_placebo)
smod_trial <- glm(S ~ X1 + X2 + trial + X1:trial + X2:trial, gaussian, df_placebo)
lrt_s <- lmtest::lrtest(smod_no_trial, smod_trial)[2,4] #test stat
lrt_s_pval <- lmtest::lrtest(smod_no_trial, smod_trial)[2,5] #p-val

lrt_y_perm <- rep(NA, n_perm); lrt_s_perm <- rep(NA, n_perm)
i <- 1
for(i in 1:n_perm){
  df_placebo_perm <- df_placebo
  df_placebo_perm$trial <- sample(df_placebo$trial)
  
  ymod_trial_perm <- glm(Y ~ X1 + X2 + trial + X1:trial + X2:trial, binomial, df_placebo_perm)
  lrt_y_perm[i] <- lmtest::lrtest(ymod_no_trial, ymod_trial_perm)[2,4] #extract LRT stat
  
  smod_trial_perm <- glm(S ~ X1 + X2 + trial + X1:trial + X2:trial, gaussian, df_placebo_perm)
  lrt_s_perm[i] <- lmtest::lrtest(smod_no_trial, smod_trial_perm)[2,4] #extract LRT stat
  i <- i + 1
}

lrt_y_perm_pval <- mean(lrt_y_perm > lrt_y)

lrt_s_perm_pval <- mean(lrt_s_perm > lrt_s)

results <- data.frame(
  lrt_y_pval, lrt_s_pval, lrt_y_perm_pval, lrt_s_perm_pval
)

output_dir <- paste(wd,"results_setting",setting,"_lrt_nt",n_t,sep="")
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

output_filename <- paste(output_dir,"/iter_",iter,".csv",sep="")
write_csv(results, output_filename)
