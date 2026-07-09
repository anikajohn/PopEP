#!/bin/bash
#SBATCH --job-name=SarsCov2_THD_clac  
#SBATCH --output=snakemake_SarsCov2.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=snakemake_SarsCov2.err        
#SBATCH --ntasks=1                  
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=2G                   
#SBATCH --partition=normal          

# Define the path to the snakemake executable
SNAKEMAKE_EXEC="/cluster/work/bewi/members/anjohn/miniconda3/envs/nextclade/bin/snakemake"

#Run per amplicon and location

# snakemake --dryrun --configfile config/haplo_SarsCov2_21_22/config.yaml > dryrun.txt

$SNAKEMAKE_EXEC  \
 --configfile config/haplo_SarsCov2_21_22/config.yaml \
 --profile profile_simple/ \
 --show-failed-logs \
 --rerun-triggers mtime \
 --rerun-incomplete 








