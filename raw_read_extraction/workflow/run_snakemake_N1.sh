#!/bin/bash
#SBATCH --job-name=N1RawReadExtraction   
#SBATCH --output=snakemake_N1.out 
#SBATCH --error=snakemake_N1.err        
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN     
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=50
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=5G                   
#SBATCH --partition=normal          

# Load the snakemake environt
# source ~/.bashrc
# conda activate lollijac

#Run per amplicon and location


# snakemake --cores 140 \
#  --configfile config/N1_24_25/config.yaml \
#  --profile profile_simple/ \
#  --show-failed-logs \
#  --until prep_samples_yaml


# snakemake all_fused_bams --cores 140 \
#  --configfile config/N1_24_25/config.yaml \
#  --profile profile_simple/ \
#  --show-failed-logs 

snakemake --cores 140 \
 --configfile config/N1_24_25/config.yaml \
 --profile profile_simple/ \
 --show-failed-logs 
