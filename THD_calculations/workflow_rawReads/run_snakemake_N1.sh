#!/bin/bash
#SBATCH --job-name=THD_calculation_N1_24_25 
#SBATCH --output=snakemake_N1.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=snakemake_N1.err        
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=50
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=5G                   
#SBATCH --partition=normal          


snakemake  \
 --configfile config/N1_24_25/config.yaml \
 --profile profile_simple/ \
 --show-failed-logs 









