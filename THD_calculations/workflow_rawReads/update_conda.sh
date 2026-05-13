#!/bin/bash
#SBATCH --job-name=conda_update
#SBATCH --time=03:30:00
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=1
#SBATCH --output=update_conda.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN

#conda before update conda 22.9.0
conda update -n base conda -y