#!/bin/bash
#SBATCH --job-name=Test_cmdl  
#SBATCH --output=cmdl_test.out         
#SBATCH --mail-user="anika.john@bsse.ethz.ch"
#SBATCH --mail-type=END,BEGIN
#SBATCH --error=cmdl_test.err     
#SBATCH --ntasks=1                  
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00             
#SBATCH --mem-per-cpu=80G                   
#SBATCH --partition=normal 


# Set variables
# output_dir="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25"
# location="over_all_locations"
# amplicon="amplicon3"
# dummy_dir=""
# reference_dataset="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/datasets/H1_26update/CY121680"
# reference="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/Lollijac/resources/H1.fasta"
# color_palette="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/resources/colour_palette/H1/colours.json"




# mkdir -p ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation
# mkdir -p ${output_dir}/${location}/${amplicon}/${dummy_dir}/subclade_abundance

# cp ${output_dir}/${location}/${amplicon}/combined_more_downsampled_amp_reads.fasta ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_more_downsampled_amp_reads.fasta
# cp ${output_dir}/${location}/${amplicon}/combined_metadata.tsv ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_metadata.tsv

# nextclade3 run \
#     --input-dataset ${reference_dataset} \
#     --input-ref ${reference} \
#     --masked-muts-weight 0 \
#     --output-all ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_results \
#     ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_more_downsampled_amp_reads.fasta

# python scripts/add_metaInfo.py \
#     -m ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_metadata.tsv \
#     -n ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_results/nextclade.tsv \
#     -l ${location} \
#     -a ${amplicon} \
#     -o ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_metadata.tsv


source /cluster/home/anjohn/miniconda3/etc/profile.d/conda.sh
conda activate THD_analysis

# # Rule 1: calculate_thd
# Rscript scripts/THD_analysis.R \
#     -q ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_more_downsampled_amp_reads.fasta \
#     -m ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_metadata.tsv \
#     -n ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_results/nextclade.tsv \
#     -r 0.00329 \
#     -s 0.15,0.2,0.25,0.3,0.35,0.4 \
#     -w 0.019,0.038,0.0833,0.15,0.23 \
#     -f ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/index_time_parameters.tsv \
#     -i ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/usedSeq_IDs.tsv \
#     -p ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/plot_index_time_parameters.png \
#     -c ${color_palette}

# # Rule 2: abundance_analysis
# Rscript scripts/abundance_analysis.R \
#     -m ${output_dir}/${location}/${amplicon}/${dummy_dir}/nextclade_metadata.tsv \
#     -i ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/usedSeq_IDs.tsv \
#     -f ${output_dir}/${location}/${amplicon}/${dummy_dir}/subclade_abundance/nextclade_metadata_processed.tsv \
#     -d ${output_dir}/${location}/${amplicon}/${dummy_dir}/subclade_abundance/plot_subclade_frequency_byDay.png \
#     -w ${output_dir}/${location}/${amplicon}/${dummy_dir}/subclade_abundance/plot_subclade_frequency_byWeek.png \
#     -c ${color_palette}
# # Rule 3: coverage_analysis
# Rscript scripts/coverage_analysis.R \
#     -c ${output_dir}/${location}/${amplicon}/combined_amp_cov_location_amplicon.tsv \
#     -f ${output_dir}/${location}/${amplicon}/amplicon_coverage/coverage_processed.tsv \
#     -b ${output_dir}/${location}/${amplicon}/amplicon_coverage/plot_covbases.png \
#     -n ${output_dir}/${location}/${amplicon}/amplicon_coverage/plot_numreads.png \
#     &> logs/${location}/${amplicon}/coverage_analysis.log


# Rscript scripts/index_visualisation.R \
#         -i ${output_dir}/${location}/${amplicon}/${dummy_dir}/THD_calculation/combined_index_time_parameters.tsv \
#         -k ${output_dir}/${location}/${amplicon}/${dummy_dir}/combined_amp_read_versions_counts.tsv \
#         -f test_index_file_processed.tsv \
#         -w test_THD_plot.png \
#         -a \
#         -c ${color_palette}




####RANDOM TEST ###########

# Rscript scripts/THD_haplo_analysis.R \
#         -q /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/combined_amp_reads_versions.fasta \
#         -k /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/combined_amp_read_versions_counts.tsv \
#         -d /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/dist_mat.tsv \
#         -m /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/combined_metadata.tsv \
#         -n /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/test_combined_nextclade_metadata.tsv \
#         -r 0.00329 \
#         -s 0.40 \
#         -w 0.019 \
#         -f /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/Test_THD_calculation/thd_0.40_0.019.tsv

# Rscript scripts/index_visualisation.R \
#             -i /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/Test_THD_calculation/thd_0.40_0.019.tsv \
#             -k /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/combined_amp_read_versions_counts.tsv \
#             -f /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/Test_THD_calculation/index_visualisation_output.tsv \
#             -w /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_H1_24_25/over_all_locations/amplicon3/Test_THD_calculation/index_visualisation_plot.png \
#             -a \
#             -c ${color_palette} \


# Rscript scripts/THD_haplo_analysis.R \
#  -q /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/raw_read_extraction/workflow/results_haplo_SarsCov2_21_22/chur/amplicon74/combined_amp_reads_versions.fasta \
#  -k /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/raw_read_extraction/workflow/results_haplo_SarsCov2_21_22/chur/amplicon74/combined_amp_read_versions_counts.tsv \
#  -d /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_SarsCov2_21_22/chur/amplicon74/dist_mat.tsv.gz \
#  -m /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_SarsCov2_21_22/chur/amplicon74/metadata.tsv \
#  -n /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_SarsCov2_21_22/chur/amplicon74/nextclade_results/nextclade.tsv \
#  -r 0.0008 \
#  -s 0.2 \
#  -w 0.23 \
#  -f /cluster/work/bewi/members/anjohn/projects/wastewater/Influenza_VaccinePred/THD_calculations/workflow_rawReads/results_haplo_SarsCov2_21_22/chur/amplicon74/THD_calculation/thd_0.2_0.23.tsv \
#  -b partiallyAliased \
#  -l clade_who


python scripts/combine_across_locations_test.py \
            --nextclade results_debug_haplo_H1_24_25/basel/amplicon4/nextclade_results/nextclade.tsv \
            results_debug_haplo_H1_24_25/chur/amplicon4/nextclade_results/nextclade.tsv \
            results_debug_haplo_H1_24_25/geneve/amplicon4/nextclade_results/nextclade.tsv \
            results_debug_haplo_H1_24_25/laupen/amplicon4/nextclade_results/nextclade.tsv \
            results_debug_haplo_H1_24_25/lugano/amplicon4/nextclade_results/nextclade.tsv \
            results_debug_haplo_H1_24_25/zurich/amplicon4/nextclade_results/nextclade.tsv  \    
            --out-nextclade results_debug_haplo_H1_24_25/over_all_locations/amplicon4/test_combined_nextclade.tsv \
            &> {log}


python scripts/add_metaInfo.py \
            -m results_debug_haplo_H1_24_25/over_all_locations/amplicon4/combined_metadata.tsv \
            -n results_debug_haplo_H1_24_25/over_all_locations/amplicon4/test_combined_nextclade.tsv \
            -l "overall_locations" \
            -a "amplicon4" \
            -o  results_debug_haplo_H1_24_25/over_all_locations/amplicon4/Test_combined_nextclade_metadata.tsv