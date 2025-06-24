# Code for "Doubly Standardized Surrogate Endpoints for SARS-CoV-2 Vaccines" 

This repository contains files for the simulation study in the manuscript "Doubly Standardized Surrogate Endpoints for SARS-CoV-2 Vaccines" by Rosin, S. P., Stivjen, F., Cross., K. A., Shook-Sa, B. E., Hudgens, M. G., & Gilbert, P. B. (2025). Simulations were conducted in R on a Slurm-managed High Performance Computing (HPC) cluster. In the below, **filenames are in bold text**. 

## Simulation settings 1-3

The shell script **sim<nolink>.sh** runs a set of simulations. The variables *setting*, *n\_t* (size of each trial), and *n\_boot* (number of bootstrap samples for estimation of $\rho_s$) are set in lines 13-15 of sim<nolink>.sh. The trial size *n\_t* used for each *setting* are given in Section 3.1 and Table S1 of the manuscript. In line 17, the array argument corresponds to a set of random seeds for the set of simulations; e.g., if "--array=1-1000" is included, 1000 simulations will be run with random seeds 1, 2, ..., 1000. In particular, the program iter.R will be run 1000 separate times with the 1000 different random seeds, with each run of iter.R correspond to one simulation iteration. The Slurm *sbatch* command is used to run sim<nolink>.sh, i.e., entering the command "sbatch sim<nolink>.sh". 

sim<nolink>.sh should be placed in the same working directory (wd) with the following five R scripts:

Conducting simulation settings 1, 2, and 3 requires the following three R scripts, which should be placed in the same working directory (wd) as sim<nolink>.sh. 
1. **iter.R** conducts a single iteration of a simulation.  **Users should edit the global variable 'wd' in line 4 to provide their working directory on their cluster.** Point estimates are saved to the file wd/setting[S]\_nt[n\_t]/iter\_[i].csv and estimated variance-covariance matrices are saved to wd/setting[S]\_nt[n\_t]/iter\_[i]_vcov.RDS, where [S], [n\_t], and [i] are numbers denoting the simulation setting, sample size in each trial, and random seed, respectively. iter.R sources the estFUN\_taudelta.R and parameter\_values.R programs.  
2. **estFUN\_taudelta.R** contains the estimating functions (equations) for variance estimation of the DSCAPs. 
3. **parameter\_values.R** contains the different values of $\theta_a$, $\gamma_a$, and $\sigma_a^2$ which control each simulation setting. The parameter values in this file for settings 1, 2, and 3 are the parameter values listed in Table S1 in the manuscript.

## Simulation setting 5 and diagnostic test from Section 2.4 of the manuscript

The diagnostic test proposed in Section 2.4 of the manuscript is evaluated with the program **lrt.R**. This program generates data analogously to iter.R, but rather than making estimates of the DSCAPs, computes likelihood ratio test and permutation test p-values for the diagnostic test. The lrt.R program is called with the shell script **lrt<nolink>.sh**, which is similar to iter.sh and allows the user to set the variables *setting*, *n_t*, and *n_perm* (number of permutations) are set. If *setting*=5, then the parameter values in **parameter_values_nonequivalent_placebo.R** are used rather than those in parameter_values.R. These values correspond to those for simulation setting 5 given in Table S2 of the manuscript, where assumption 6 of conditional exchangeability of trial did not hold for the placebo group.

## Computing results

In **truevalue.R**, the true value of each DSCAP was determined empirically based on the average potential outcomes from $2 \times 10^7$ individuals randomly sampled from the target population. This script is only run once for each simulation setting, but is called from the shell script **truevalue<nolink>.sh**. Set the two parameters *setting* and *n_t* in truevalue<nolink>.sh; *n_t* should be set to $2 \times 10^7$. 

For a given simulation setting, once the true value of the DSCAPs have been determined with truevalue.R, and the simulations have been conducted using iter.R, sim<nolink>.sh, and lrt.R, the results  (bias, confidence interval coverage, etc.) can be computed from **results.R**. The simulation setting and corresponding *n_t* and number of simulations *n_sim* are set in the first three lines of results.R. 

