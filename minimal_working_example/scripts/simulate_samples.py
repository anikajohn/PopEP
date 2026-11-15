"""Simulate amplicon read pairs per sample with art_illumina and align them to the reference.

For every location x sampling date:
  1. expected version frequencies from the scenario dynamics
  2. per amplicon: multinomial draw of read_pairs_per_amplicon fragments over the versions
  3. art_illumina (amplicon mode, paired-end) on the fragments; read names keep the version
  4. optional edge case: lower base quality at a variable site in a fraction of reads
  5. bwa mem -> sorted, indexed BAM in the V-pipe layout PopEP expects:
       <outdir>/samples/<A1>_<code>_<YYYY>_<MM>_<DD>/<YYYYMMDD>_SIM/alignments/REF_aln_trim.bam

Truth tables go to --truthdir: expected_frequencies.tsv, fragment_counts.tsv,
low_quality_reads.tsv, samples.tsv.
"""
import argparse
import datetime as dt
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from Bio import SeqIO

from common import (batch_name, expected_frequencies, load_scenario, read_inserts,
                    read_reference, resolve, sample_name, sampling_dates, _date)


def run(cmd, **kw):
    subprocess.run(cmd, check=True, **kw)


def read_fastq(path):
    with open(path) as f:
        while True:
            head = f.readline()
            if not head:
                return
            seq, plus, qual = f.readline(), f.readline(), f.readline()
            yield [head.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")]


def write_fastq(records, handle):
    for r in records:
        handle.write("\n".join(r) + "\n")


def fragment_id(read_header):
    """ART names reads '<fragment id>-<n>/<mate>'; return the fragment id."""
    return read_header[1:].split()[0].rsplit("/", 1)[0].rsplit("-", 1)[0]


def low_quality_rules(cfg, location, date):
    rules = []
    for r in cfg.get("low_quality_sites", []) or []:
        if location in r["locations"] and _date(r["from"]) <= date <= _date(r["to"]):
            rules.append(r)
    return rules


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--haplotypes", required=True, help="haplotypes.fasta from build_haplotypes.py")
    ap.add_argument("--outdir", required=True, help="where the sample folders are written")
    ap.add_argument("--truthdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--keep-fastq", action="store_true")
    args = ap.parse_args()

    cfg = load_scenario(args.config)
    seq_cfg = cfg["sequencing"]
    ref_id, _ = read_reference(cfg)
    ref_path = resolve(cfg["reference"])
    inserts = read_inserts(cfg)
    haps = {r.id: str(r.seq) for r in SeqIO.parse(resolve(args.haplotypes), "fasta")}
    names = list(cfg["haplotypes"])
    rng = np.random.default_rng(cfg["seed"])
    pad, read_len, n_pairs = seq_cfg["fragment_padding"], seq_cfg["read_length"], seq_cfg["read_pairs_per_amplicon"]

    outdir, truthdir = resolve(args.outdir), resolve(args.truthdir)
    sample_root = outdir / "samples"
    if sample_root.exists():
        shutil.rmtree(sample_root)
    sample_root.mkdir(parents=True)
    truthdir.mkdir(parents=True, exist_ok=True)

    # bwa index lives next to a private copy of the reference, so resources/ stays clean
    idx_dir = outdir / "bwa_index"
    idx_dir.mkdir(exist_ok=True)
    idx_ref = idx_dir / ref_path.name
    shutil.copy(ref_path, idx_ref)
    run(["bwa", "index", str(idx_ref)], stderr=subprocess.DEVNULL)

    freq_out = open(truthdir / "expected_frequencies.tsv", "w")
    count_out = open(truthdir / "fragment_counts.tsv", "w")
    lq_out = open(truthdir / "low_quality_reads.tsv", "w")
    samples_out = open(truthdir / "samples.tsv", "w")
    freq_out.write("location\tdate\thaplotype\texpected_frequency\n")
    count_out.write("sample\tlocation\tdate\tamplicon\thaplotype\tn_fragments\n")
    lq_out.write("sample\tamplicon\tfragment\tposition\tquality\n")
    samples_out.write("sample\tbatch\tlocation\tlocation_code\tlocation_name\tdate\tbam\n")

    art_seed = int(cfg["seed"])
    for location, loc in cfg["locations"].items():
        for date in sampling_dates(cfg):
            sample, batch = sample_name(cfg, location, date), batch_name(date)
            freqs = expected_frequencies(cfg, location, date)
            for h in names:
                freq_out.write(f"{location}\t{date}\t{h}\t{freqs[h]:.6f}\n")
            p = np.array([freqs[h] for h in names])

            with tempfile.TemporaryDirectory() as tmp:
                tmp = Path(tmp)
                r1_all, r2_all = open(tmp / "R1.fq", "w"), open(tmp / "R2.fq", "w")
                for amp, (s, e) in inserts.items():
                    fs, fe = s - pad, e + pad
                    counts = rng.multinomial(n_pairs, p)
                    frag_fa = tmp / f"{amp}.fa"
                    with open(frag_fa, "w") as f:
                        for h, c in zip(names, counts):
                            count_out.write(f"{sample}\t{location}\t{date}\t{amp}\t{h}\t{c}\n")
                            for i in range(c):
                                f.write(f">{h}__{amp}__{i}\n{haps[h][fs:fe]}\n")
                    art_seed += 1
                    cmd = ["art_illumina", "-ss", seq_cfg["art_profile"], "-amp", "-p", "-na", "-q",
                           "-l", str(read_len), "-f", "1", "-rs", str(art_seed),
                           "-i", str(frag_fa), "-o", str(tmp / f"{amp}_")]
                    if seq_cfg.get("no_indel_errors", True):
                        cmd += ["-ir", "0", "-ir2", "0", "-dr", "0", "-dr2", "0"]
                    if seq_cfg.get("min_base_quality", 0):
                        cmd += ["-qL", str(seq_cfg["min_base_quality"])]
                    # stderr only carries "will not output any ALN or SAM file" (intended: we align ourselves)
                    run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    r1 = list(read_fastq(tmp / f"{amp}_1.fq"))
                    r2 = list(read_fastq(tmp / f"{amp}_2.fq"))

                    # Edge case: low base quality at a variable site in a fraction of fragments
                    frag_len = fe - fs
                    for rule in low_quality_rules(cfg, location, date):
                        pos0 = int(rule["position"]) - 1
                        if not (s <= pos0 < e):
                            continue
                        off = pos0 - fs
                        qchar = chr(int(rule["quality"]) + 33)
                        hit = rng.random(len(r1)) < float(rule["fraction_of_reads"])
                        for k in np.flatnonzero(hit):
                            for rec, o in ((r1[k], off), (r2[k], frag_len - 1 - off)):
                                if 0 <= o < len(rec[3]):
                                    rec[3] = rec[3][:o] + qchar + rec[3][o + 1:]
                            lq_out.write(f"{sample}\t{amp}\t{fragment_id(r1[k][0])}\t{rule['position']}\t{rule['quality']}\n")
                    write_fastq(r1, r1_all)
                    write_fastq(r2, r2_all)
                r1_all.close()
                r2_all.close()

                aln_dir = sample_root / sample / batch / "alignments"
                aln_dir.mkdir(parents=True)
                bam = aln_dir / "REF_aln_trim.bam"
                bwa = subprocess.Popen(["bwa", "mem", "-v", "1", "-t", str(args.threads), str(idx_ref),
                                        str(tmp / "R1.fq"), str(tmp / "R2.fq")],
                                       stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
                sort = subprocess.Popen(["samtools", "sort", "--no-PG", "-o", str(bam), "-"],
                                        stdin=subprocess.PIPE, text=True)
                for line in bwa.stdout:
                    if not line.startswith("@PG"):
                        sort.stdin.write(line)
                sort.stdin.close()
                if bwa.wait() or sort.wait():
                    raise RuntimeError(f"alignment failed for {sample}")
                run(["samtools", "index", str(bam)])
                if args.keep_fastq:
                    shutil.copy(tmp / "R1.fq", aln_dir.parent / "R1.fq")
                    shutil.copy(tmp / "R2.fq", aln_dir.parent / "R2.fq")

            samples_out.write(f"{sample}\t{batch}\t{location}\t{loc['code']}\t{loc['name']}\t{date}\t"
                              f"{bam.relative_to(outdir)}\n")
            print(f"{sample}: done")

    for f in (freq_out, count_out, lq_out, samples_out):
        f.close()


if __name__ == "__main__":
    main()
