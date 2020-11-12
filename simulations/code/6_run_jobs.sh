#!/bin/bash

#PBS -N sim_sfpca
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=10gb
#PBS -V
#PBS -j oe
#PBS -d .
#PBS -t 100
#PBS -o messages_outputs/
#PBS -e messages_errors/

set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER/sfpca_2020
[ ! -d $TMPDIR ] && mkdir $TMPDIR
export TMPDIR=$TMPDIR/simulations
[ ! -d $TMPDIR ] && mkdir $TMPDIR

# load module
module load R_3.6.3-9.3.0

# simulate data (t: number of each simulation scenario, e.g. 1-16%10)
#Rscript 1_generate_simulated_data.R $TMPDIR

# apply BayesTime (t: number of simulation index: #replicates * scenarios i.e. 1k * 16 # barnacle maximum 1k each)
#Rscript 2_bayesTime_simulations.R $TMPDIR
#Rscript 5_model_selection_diagnostics.R $TMPDIR

# summarize results (t: number of replicates for each scenario, i.e. 1k)
Rscript 3_summary_sim_results_quantiles.R $TMPDIR


