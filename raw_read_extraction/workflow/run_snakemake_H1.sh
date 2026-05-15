#!/bin/bash
#SBATCH --job-name=H1RawReadExtraction   
#SBATCH --output=snakemake.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=snakemake.err        
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=50
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=5G                   
#SBATCH --partition=normal          

# Load the snakemake environt
# source ~/.bashrc
# conda activate lollijac

#Run per amplicon and location

# snakemake --dryrun --configfile config/H1_Geneva_05Jan25/config.yaml > dryrun.txt


snakemake --cores 140 \
 --configfile config/H1_25_26/config.yaml \
 --profile profile_simple/ \
 --show-failed-logs 



