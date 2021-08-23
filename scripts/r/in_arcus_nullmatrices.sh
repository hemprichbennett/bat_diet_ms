#!/bin/bash

#SBATCH --job-name=null_models # the name for the cluster scheduler
#SBATCH --time=4:00:00 # Maximum allowed runtime per iteration
#SBATCH --array=1-1000 # the number of iterations
#SBATCH --mem-per-cpu=8G   # memory per cpu-core
#SBATCH --output=arcus_outputs/null_models%A_%a.out # the name of the output files
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.hemprich-bennett@zoo.ox.ac.uk

module load python/anaconda3/2019.03

source activate $DATA/conda_envs/bat_diet

Rscript arcus_null_matrix_gen.R ${SLURM_ARRAY_TASK_ID}
