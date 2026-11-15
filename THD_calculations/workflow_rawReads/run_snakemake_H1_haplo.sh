#!/bin/bash
#SBATCH --job-name=H1_THD_clac  
#SBATCH --output=snakemake_H1_haplo.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=snakemake_H1_haplo.err        
#SBATCH --ntasks=1                  
#SBATCH --time=72:00:00             
#SBATCH --mem-per-cpu=2G                   
#SBATCH --partition=normal          

# Define the path to the snakemake executable
SNAKEMAKE_EXEC="/cluster/work/bewi/members/anjohn/miniconda3/envs/nextclade/bin/snakemake"

#Run per amplicon and location

# snakemake --dryrun --configfile config/haplo_H1_24_25/config.yaml  --rerun-triggers mtime --rerun-incomplete > dryrun.txt

$SNAKEMAKE_EXEC  \
 --configfile config/haplo_H1_24_25/config.yaml \
 --profile profile_simple/ \
 --show-failed-logs \
 --rerun-triggers mtime \
 --rerun-incomplete 








