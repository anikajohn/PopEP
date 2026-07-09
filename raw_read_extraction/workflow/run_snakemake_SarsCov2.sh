#!/bin/bash
#SBATCH --job-name=SarsCov2HaploReads  
#SBATCH --output=SarsCov2HaploReads.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=SarsCov2HaploReads.err        
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=50
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=5G                   
#SBATCH --partition=normal          

# Load the snakemake environt
# source ~/.bashrc
# conda activate lollijac

#Run per amplicon and location


snakemake --cores 140 \
 --configfile config/SarsCov2_21_22/config.yaml \
 --profile profile_simple/ \
 --show-failed-logs \
 --rerun-incomplete \
 --rerun-triggers mtime



