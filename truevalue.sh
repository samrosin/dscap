#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 01:00:00
#SBATCH --mem=32g
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=srosin@bsc.gwu.edu
module load r/4.1.0

# define variables
setting=1
n_t=2e+07

Rscript truevalue.R $setting $n_t 
