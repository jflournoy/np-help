#!/bin/bash
#Slurm submission options
#LOOK AT THESE AND EDIT TO OVERRIDE FOR YOUR JOB
#SBATCH --job-name buildimg
#SBATCH -p ncf_holy
#SBATCH --mem 40000
#SBATCH --cpus-per-task 64
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --time 2-00:00
#SBATCH --mail-type=END
#SBATCH -o %x_%A.out

module load gcc/7.1.0-fasrc01
module load R/3.5.1-fasrc01

PREFIX=$1

export OMP_NUM_THREADS=64
Rscript ~/otherhome/code/np-help/build_image_from_rds.R --prefix "${PREFIX}" --ncores 64 $2