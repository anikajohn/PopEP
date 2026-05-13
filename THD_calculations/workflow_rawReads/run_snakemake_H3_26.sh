#!/bin/bash
#SBATCH --job-name=H3THD_calculation 
#SBATCH --output=snakemake_H3_26.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=snakemake_H3_26.err        
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=50
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=5G                   
#SBATCH --partition=normal          


snakemake  \
 --configfile config/H3_25_26/config.yaml \
 --profile profile_simple/ \
 --show-failed-logs 









