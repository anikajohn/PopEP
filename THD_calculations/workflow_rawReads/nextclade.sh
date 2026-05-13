#!/bin/bash
#SBATCH --job-name=Test_cmdl  
#SBATCH --output=cmdl_test.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=cmdl_test.err     
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=50
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=10G                   
#SBATCH --partition=normal 

# Set variables
output_dir="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/scripts/clinical_data/H1Gisaid_2425_2526/nextclade/"
reference_dataset="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/datasets/H1_26update/CY121680"
reference="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/Lollijac/resources/H1.fasta"


nextclade3 run \
    --input-dataset ${reference_dataset} \
    --input-ref ${reference} \
    --masked-muts-weight 0 \
    --output-all ${output_dir}/nextclade_results \
    ${output_dir}/sequence.fasta
