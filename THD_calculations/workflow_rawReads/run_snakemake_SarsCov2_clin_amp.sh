#!/bin/bash
#SBATCH --job-name=AmpSarsCov2_THD_clac  
#SBATCH --output=snakemake_SarsCov2_amp.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=snakemake_SarsCov_amp2.err        
#SBATCH --ntasks=1                  
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=2G                   
#SBATCH --partition=normal          

# Define the path to the snakemake executable
SNAKEMAKE_EXEC="/cluster/work/bewi/members/anjohn/miniconda3/envs/nextclade/bin/snakemake"

#Run per amplicon and location

# snakemake --dryrun --configfile config/clin_SarsCov2_21_22/config.yaml > dryrun.txt

$SNAKEMAKE_EXEC  \
 --configfile config/clin_SarsCov2_21_22/config.yaml \
 --profile profile_simple/ \
 --show-failed-logs \
 --rerun-triggers mtime \
 --rerun-incomplete 








