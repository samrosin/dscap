#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 05:00:00
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=srosin@bsc.gwu.edu
module load r/4.1.0 

#define variables
setting=1
n_t=3000
n_perm=10000

sbatch --output=/dev/null --error=/dev/null --time=02:00:00 --mem=1g --array=10001-11020 --job-name=lrt --wait R CMD BATCH --no-save --no-restore "--args $setting $n_t $n_perm" lrt.R lrt.Rout
