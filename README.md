# dscap
Doubly standardized causal association parameters for correlate of protection evaluation.

These programs are for running simulations corresponding to the manuscript "Doubly Standardized Surrogate Endpoints for SARS-CoV-2 Vaccines". The simulations are intended to be run on a Slurm-managed High Performance Computing (HPC) cluster, such as UNC's Longleaf cluster. 

Each simulation iteration is run in iter.R. The only line of code that should need to be changed is line 10, which sets the working directory. The iter.R program sources estFUN_taudelta.R, which contains code for the estimating equations for the regression slope and Pearson correlation DSCAPs, and parameter_values.R. 

The set of simulations is run from the sim.sh shell script. There are three editable parameters in sim.sh, in lines 13-15: setting, which corresponds to simulation setting; n_t, which corresponds to the number of participants in each trial; and n_boot, the number of bootstrap samples, which can be ignored for the purposes of the regression slope and Pearson correlation DSCAPs. The number of simulations is controlled by the "--array=[#]-[#]" parameter in the sbatch statement on line 17, where each [#] can be replaced with a range of numbers, e.g., a range such as 1-1000 would correspond to 1000 simulations (run with random seeds 1 through 1000). The time and memory allocated for each job can also be edited in the sim.sh file. The command "sbatch sim.sh" is run from a cluster, e.g., the Longleaf cluster, in order to perform the simulations. 

The results.R program can be used to analyze the bias and confidence interval coverage of the estimators across a set of simulations. 
