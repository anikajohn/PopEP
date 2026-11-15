# Minimal example: simulated H1 wastewater amplicon data

Simulated input for the PopEP `raw_read_extraction` and `THD_calculations` workflows,
with known version frequencies so pipeline results can be checked against ground truth.

## Run

```bash
conda env create -f environment.yaml        # once
conda activate popep_sim
bash run_simulation.sh 8                    # 8 threads, ~1 min
```

Everything is driven by `config/scenario.yaml` and is deterministic (fixed `seed`).

## Scenario

- Reference: CY121680.1 (`resources/reference/H1.fasta`), amplicons 2–6 (`resources/amplicons/`).
- Locations: Zürich (10), Lugano (05); weekly samples 2024-11-01 to 2025-04-01 (22 per location).
- Versions (built from `C193_example.fasta` / `D31_example.fasta`, projected onto the reference):

  | version | definition | dynamics |
  |---|---|---|
  | C.1.9.3_v1 | example sequence | weight 0.5, constant |
  | C.1.9.3_v2 | v1 + C400T (amp2) | weight 0.3, constant |
  | C.1.9.3_v3 | v2 + G1100A (amp4, amp5) | weight 0.2, constant |
  | D.3.1_v1 | example sequence | from D.3.1 start: 0.02 · e^(0.08·days) |
  | D.3.1_v2 | v1 + C1540T (amp6) | from D.3.1 start: 0.02 · e^(0.05·days) |

  Frequencies are the weights renormalised, so all C.1.9.3 versions decline at the same speed.
  D.3.1 starts 2025-02-01 in Zürich and 2025-02-14 in Lugano.
- Reads: per sample and amplicon, 200 fragments drawn multinomially from the expected
  frequencies; fragment = insert ± 5 bp; `art_illumina` amplicon mode, 2×250 bp, MSv3 profile,
  no indel errors, minimum base quality 25 (so the only low-quality bases are the edge case).
- Edge case: in Zürich samples from 2024-12-01 to 2025-01-15, 50% of reads covering position 400
  get base quality 2 there. PopEP's `scan_bam.py` treats Q ≤ `phred_score_threshold` (5 in the current H1 configs) as not covered, so
  `simulate_reads.py` writes `N` at 400 for those reads.

## Outputs

| path | content |
|---|---|
| `data/samples/A1_<code>_<YYYY>_<MM>_<DD>/<YYYYMMDD>_SIM/alignments/REF_aln_trim.bam` | input for PopEP (`path_to_data: data/samples`) |
| `results/truth/haplotypes.fasta` | version sequences, reference coordinates |
| `results/truth/haplotype_mutations.tsv`, `pairwise_distances.tsv` | differences to reference / between versions |
| `results/truth/expected_frequencies.tsv` | expected frequency per location, date, version |
| `results/truth/fragment_counts.tsv` | simulated fragments per sample, amplicon, version |
| `results/truth/low_quality_reads.tsv` | reads that received the low-quality base |
| `results/check/observed_abundance.tsv` | read counts per version recovered from the BAMs vs expected |
| `results/check/abundance_over_time.png`, `clade_over_time.png` | expected (lines) vs observed (points) |

Read names carry the true version: `<version>__<amplicon>__<i>-<n>/<mate>`.

## Running PopEP raw_read_extraction on it (`pipeline_test/`)

```bash
conda env create -f pipeline_test/environment_raw_read_extraction.yaml   # once
conda activate popep_raw
bash pipeline_test/run_raw_read_extraction.sh /path/to/PopEP 8          # writes results/pipeline/
conda activate popep_sim
python pipeline_test/compare_pipeline_output.py --config config/scenario.yaml \
    --truthdir results/truth --pipeline results/pipeline/raw_read_extraction --outdir results/pipeline_check
```

The comparison merges versions that are identical within an amplicon, assigns every read version
to the closest group and writes `pipeline_vs_truth.tsv`, `read_versions_assigned.tsv` and
`low_quality_N_check.tsv`.

### Both workflows (raw_read_extraction → THD_calculations)

THD_calculations calls conda envs by name (`conda run -n nextclade` / `-n THD_analysis`). Test copies
are built into a separate envs folder, which the run script puts first on `CONDA_ENVS_PATH`, so
your own `nextclade` env is neither used nor changed:

```bash
conda env create -p ~/miniconda3/envs_popep_test/nextclade    -f pipeline_test/environment_nextclade.yaml
conda env create -p ~/miniconda3/envs_popep_test/THD_analysis -f pipeline_test/environment_THD_analysis.yaml
conda activate popep_raw
bash pipeline_test/run_both_workflows.sh /path/to/PopEP 8     # ~12 min; results/pipeline/THD_calculations/
```
