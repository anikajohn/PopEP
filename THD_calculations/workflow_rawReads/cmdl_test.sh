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
output_dir="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_H1_25_26"
location="over_all_locations"
amplicon="amplicon5"
dummy_dir="updated_ref_dataset"
reference_dataset="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/datasets/H1_26update/CY121680"
reference="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/Lollijac/resources/H1.fasta"
color_palette="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/colour_palette/H1/colours.json"

mkdir -p ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation
mkdir -p ${output_dir}/${location}/${amplicon}/${dummy_dir}/subclade_abundance

cp ${output_dir}/${location}/${amplicon}/combined_more_downsampled_amp_reads.fasta ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_more_downsampled_amp_reads.fasta
cp ${output_dir}/${location}/${amplicon}/combined_metadata.tsv ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_metadata.tsv

nextclade3 run \
    --input-dataset ${reference_dataset} \
    --input-ref ${reference} \
    --masked-muts-weight 0 \
    --output-all ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_results \
    ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_more_downsampled_amp_reads.fasta

python scripts/add_metaInfo.py \
    -m ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_metadata.tsv \
    -n ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_results/nextclade.tsv \
    -l ${location} \
    -a ${amplicon} \
    -o ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_metadata.tsv


source /cluster/home/anjohn/miniconda3/etc/profile.d/conda.sh
conda activate THD_analysis

# Rule 1: calculate_thd
Rscript scripts/THD_analysis.R \
    -q ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_more_downsampled_amp_reads.fasta \
    -m ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_metadata.tsv \
    -n ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_results/nextclade.tsv \
    -r 0.00329 \
    -s 0.15,0.2,0.25,0.3,0.35,0.4 \
    -w 0.019,0.038,0.0833,0.15,0.23 \
    -f ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/index_time_parameters.tsv \
    -i ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/usedSeq_IDs.tsv \
    -p ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/plot_index_time_parameters.png \
    -c ${color_palette}

# Rule 2: abundance_analysis
Rscript scripts/abundance_analysis.R \
    -m ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_metadata.tsv \
    -i ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/usedSeq_IDs.tsv \
    -f ${output_dir}/${location}/${amplicon}/${dummy_dir}/subclade_abundance/nextclade_metadata_processed.tsv \
    -d ${output_dir}/${location}/${amplicon}/${dummy_dir}/subclade_abundance/plot_subclade_frequency_byDay.png \
    -w ${output_dir}/${location}/${amplicon}/${dummy_dir}/subclade_abundance/plot_subclade_frequency_byWeek.png \
    -c ${color_palette}
# # Rule 3: coverage_analysis
# Rscript scripts/coverage_analysis.R \
#     -c ${output_dir}/${location}/${amplicon}/combined_amp_cov_location_amplicon.tsv \
#     -f ${output_dir}/${location}/${amplicon}/amplicon_coverage/coverage_processed.tsv \
#     -b ${output_dir}/${location}/${amplicon}/amplicon_coverage/plot_covbases.png \
#     -n ${output_dir}/${location}/${amplicon}/amplicon_coverage/plot_numreads.png \
#     &> logs/${location}/${amplicon}/coverage_analysis.log