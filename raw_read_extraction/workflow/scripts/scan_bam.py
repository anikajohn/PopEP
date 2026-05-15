import pandas as pd
import pysam
import json
"""
Analyzes the alignment file (BAM) of fused amplicons.
Input: alignment file, reference, BED insert file, amplicons dictionary.
Parameters: PHRED score threshold.
Output:
dictionary:
keys - numbers of amplicon regions ( 1 to 98 ),
values - dictionary {
            tuple of combination: {
                "amplicon": amp_number,
                "covered": covered_amp_binary,
                "mutated": mutated_amp_binary,
                "counts": 54        <-- raw read counts (not normalised)
                }
}
"""
print("Inputs:", snakemake.input)
print("BAM path:", snakemake.input.fused_bam)
alnfname = snakemake.input.fused_bam  # bam file with FUSED reads (i.e. AMPLICONS)
amplicons_dict_file = snakemake.input.amplicons_dict
output = snakemake.output.output
reference_name = snakemake.params.reference

with open(amplicons_dict_file) as f:
    amplicons_dict = json.load(f)

rq_chr = reference_name
phred_threshold_param = snakemake.params.phred


def analyse_amplicon(amplicon, amp_dict, phred_thres=phred_threshold_param):
    """
    Test if mutations (listed in amp_dict) are present on the specific pysam amplicon.
    Returns a tuple of 2 lists with:
            covered_sites:  list of sites covered (with quality > phred_thres) on the amplicon (no matter which nt)
            found_mutations: list of those positions which have the mutated variant
    """
    # WARNING pysam is 0-based! (see here: https://pysam.readthedocs.io/en/latest/faq.html#pysam-coordinates-are-wrong)
    # 1. check which positions are covered on that amplicon
    aligned = amplicon.get_aligned_pairs(matches_only=True)
    aligned_ref_to_query_pos = {r: q for q, r in aligned}  # 0-based

    covered_sites = []
    for pos in amp_dict.keys():  # 1-based
        pysam_pos = int(pos) - 1
        if pysam_pos in aligned_ref_to_query_pos.keys():
            if amplicon.query_qualities[aligned_ref_to_query_pos[pysam_pos]] > phred_thres:
                covered_sites.append(pos)

    if not len(covered_sites):
        return None, None  # sites quality is too low or position not in this read

    found_mutations = []
    for pos in covered_sites:  # covered_sites is 1-based
        pysam_pos = int(pos) - 1
        muts_at_pos = amp_dict.get(pos, [])
        if not muts_at_pos:
            continue
        hit_list = [pos if amplicon.query_sequence[aligned_ref_to_query_pos[pysam_pos]] == mut else 0
                    for mut in muts_at_pos]

        if any(hit_list):
            if len(hit_list) == 1:
                found_mutations.append(hit_list[0])
            else:
                found_mutations.append(hit_list)

    if len(found_mutations):
        return covered_sites, found_mutations
    else:
        return covered_sites, None


def pos_in_found(pos, t_mutated):
    for x in t_mutated:
        if isinstance(x, list):
            if pos in x:
                return True
        elif x == pos:
            return True
    return False


def extract_inner_list(pos, t_mutated):
    for x in t_mutated:
        if isinstance(x, list):
            if pos in x:
                return x


def freeze(obj):
    """Recursively convert lists -> tuples so object becomes hashable."""
    if isinstance(obj, list):
        return tuple(freeze(x) for x in obj)
    return obj


def scanamplicon(amp_number, amps_from_amp_region, amp_dict):
    """
    Input:
    amp_number - which INSERT REGIONS (1 to 98) is being ANALYSED
    amps_from_amp_region - amplicons that are within the region of the amplicon
    amp_dict - mutations dictionary {position: mutated_nt, ...}
    """
    amplicon_info = {}

    amplicons = {}
    for amplicon in amps_from_amp_region:
        name = str(amplicon.query_name)
        amplicons[name] = amplicon

    for amp in amplicons.values():
        amp_dict_keys = list(amp_dict.keys())
        covered_amp_binary = [0] * len(amp_dict_keys)
        mutated_amp_binary = [0] * len(amp_dict_keys)

        t_covered, t_mutated = analyse_amplicon(amp, amp_dict)

        for i, pos in enumerate(amp_dict_keys):
            if t_covered is not None and pos in t_covered:
                covered_amp_binary[i] = 1

            muts_at_pos = amp_dict.get(pos, [])
            if t_mutated is not None and pos_in_found(pos, t_mutated):
                if len(muts_at_pos) > 1:
                    inner_list = extract_inner_list(pos, t_mutated)
                    mutated_amp_binary[i] = [
                        1 if int(p) > 0 else 0
                        for p in inner_list]
                else:
                    mutated_amp_binary[i] = 1
            elif t_mutated is not None and not pos_in_found(pos, t_mutated):
                if len(muts_at_pos) > 1:
                    mutated_amp_binary[i] = [0 for p in muts_at_pos]
                else:
                    mutated_amp_binary[i] = 0

        covered_frozen = freeze(covered_amp_binary)
        mutated_frozen = freeze(mutated_amp_binary)

        combination_key = (covered_frozen, mutated_frozen)

        if combination_key in amplicon_info:
            amplicon_info[combination_key]["counts"] += 1
        else:
            amplicon_info[combination_key] = {
                "amplicon": amp_number,
                "covered": covered_amp_binary,
                "mutated": mutated_amp_binary,
                "counts": 1  # raw count, not normalised
            }

    return amplicon_info


# Analyze the whole bam file
# Note: amp_fused_sorted.bam is already filtered to this amplicon by extract_amplicon_reads,
# so no coordinate-based fetch or position filtering is needed.
amp_results = {}

with pysam.AlignmentFile(alnfname, "rb") as alnfile:
    if rq_chr is None:
        raise ValueError("Reference not provided")

    for key, values in amplicons_dict.items():
        if values is None:
            continue

        all_reads = list(alnfile.fetch(reference=rq_chr))

        single_amp_results = scanamplicon(
            amp_number=key,
            amps_from_amp_region=all_reads,
            amp_dict=amplicons_dict[key]
        )

        # Store raw counts directly -- no normalisation by total reads
        amp_results[key] = single_amp_results

amp_results_df = pd.DataFrame([
    {"key": k, "value": v} for k, v in amp_results.items()
])
amp_results_df.to_csv(output, header=None, index=False)
