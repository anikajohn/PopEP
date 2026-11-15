#!/bin/bash
#SBATCH --job-name=AmpCov_SarsCov2  
#SBATCH --output=AmpCov_SarsCov2.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=AmpCov_SarsCov2.err        
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00             
#SBATCH --mem-per-cpu=5G                   
#SBATCH --partition=normal    

# python cov_amplicon_insert.py \
#   -b v3=/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/raw_read_extraction/workflow/resources/SarsCov2_21_22/primers/v3/nCoV-2019.insert.bed \
#   -b v4=/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/raw_read_extraction/workflow/resources/SarsCov2_21_22/primers/v4/SARS-CoV-2.insert.bed \
#   -s target_samples.tsv \
#   -p /cluster/project/pangolin/processes/sars_cov_2/vpipe_output \
#   -o AmpCov_SarsCov2_results 

python cov_amplicon_insert.py \
  -b dummy=/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/raw_read_extraction/workflow/resources/H1_24_25/dummy_primers/H1dummy.insert.bed \
  -s H1_samples.tsv \
  -p /cluster/project/pangolin/old_setup/not_necessary/influenza_pipeline/IA_H1/working/results \
  -o AmpCov_H1_results 
