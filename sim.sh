#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 72:00:00
#SBATCH --mem=5g
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=srosin@live.unc.edu
module load r/4.1.0 

#define variables
setting=3
n_t=8000
n_boot=2 # number of bootstrap samples -- irrelevant for sim_beta

sbatch --output=/dev/null --error=/dev/null --time=48:00:00 --mem=5g --array=60001-60005 --job-name=iter --wait R CMD BATCH --no-save --no-restore "--args $setting $n_t $n_boot" iter.R iter.Rout
