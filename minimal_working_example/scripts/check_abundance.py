"""Check the simulated BAMs against the scenario.

Reads every sample BAM, takes read 1 of each mapped primary pair (MAPQ >= 9, as in PopEP),
recovers its version from the read name and compares observed with expected frequencies.
Also checks the low-quality edge case and the per-read error rate against the true version.

Outputs (in --outdir):
  observed_abundance.tsv      location, date, amplicon, haplotype, reads, observed vs expected
  low_quality_check.tsv       per sample: reads covering the site, reads with Q <= 20 there
  abundance_over_time.png     expected (lines) vs observed (points) per version and location
  clade_over_time.png         same, summed per clade
"""
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import pysam
from Bio import SeqIO

from common import _date, clade_of, load_scenario, read_inserts, resolve

MIN_MAPQ = 9          # PopEP minimal_read_mapping_quality
PHRED_THRESHOLD = 5   # overwritten by scenario pipeline_phred_threshold (scan_bam.py: covered if quality > it)

# Reference categorical palette, fixed order (dataviz skill, light mode); marker = secondary encoding
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
MARKERS = ["o", "s", "^", "D", "v"]
INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"


def count_reads(bam_path, haps, lq_positions):
    counts = defaultdict(int)
    lq = defaultdict(lambda: [0, 0])  # position -> [covering, low quality]
    mism, bases = 0, 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped or r.is_secondary or r.is_supplementary or not r.is_read1 \
                    or r.mate_is_unmapped or r.mapping_quality < MIN_MAPQ:
                continue
            hap, amp, _ = r.query_name.split("__")
            counts[(amp, hap)] += 1
            seq, qual, truth = r.query_sequence, r.query_qualities, haps[hap]
            for q, ref in r.get_aligned_pairs(matches_only=True):
                bases += 1
                mism += seq[q] != truth[ref]
                if ref + 1 in lq_positions:
                    lq[ref + 1][0] += 1
                    lq[ref + 1][1] += qual[q] <= PHRED_THRESHOLD
    return counts, lq, mism, bases


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK_2)
    ax.tick_params(colors=INK_2, labelsize=9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %y"))
    ax.set_ylim(-0.02, 1.02)


def plot(df, exp, cfg, group, path, title):
    names = list(dict.fromkeys(exp[group]))
    locs = list(cfg["locations"])
    fig, axes = plt.subplots(1, len(locs), figsize=(6 * len(locs), 4.2), sharey=True, facecolor=SURFACE)
    for ax, loc in zip(axes, locs):
        style(ax)
        start = pd.Timestamp(_date(cfg["locations"][loc]["d31_start"]))
        ax.axvline(start, color=INK_2, linestyle="--", linewidth=1)
        ax.text(start, 1.0, " D.3.1 start", color=INK_2, fontsize=8, va="top")
        for rule in cfg.get("low_quality_sites", []) or []:
            if loc in rule["locations"]:
                ax.axvspan(pd.Timestamp(_date(rule["from"])), pd.Timestamp(_date(rule["to"])),
                           color=GRID, alpha=0.6, linewidth=0)
                ax.text(pd.Timestamp(_date(rule["from"])), 0.93, f" low-Q at {rule['position']}",
                        color=INK_2, fontsize=8, va="top")
        for i, n in enumerate(names):
            e = exp[(exp.location == loc) & (exp[group] == n)].groupby("date").expected_frequency.sum()
            o = df[(df.location == loc) & (df[group] == n)].groupby(["date", "amplicon"]).observed_frequency.sum()
            o = o.groupby("date").agg(["mean", "min", "max"])
            ax.plot(e.index, e.values, color=SERIES[i], linewidth=2, label=n)
            ax.errorbar(o.index, o["mean"], yerr=[o["mean"] - o["min"], o["max"] - o["mean"]], fmt=MARKERS[i],
                        color=SERIES[i], markersize=6, markeredgecolor=SURFACE, markeredgewidth=1,
                        elinewidth=1, capsize=0)
        ax.set_title(f"{cfg['locations'][loc]['name']}", color=INK, fontsize=11, loc="left")
    axes[0].set_ylabel("Frequency", color=INK_2)
    handles = [plt.Line2D([], [], color=SERIES[i], marker=MARKERS[i], linewidth=2, markersize=6) for i in range(len(names))]
    fig.legend(handles, names, loc="upper center", ncol=len(names), frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0.99), labelcolor=INK)
    fig.suptitle(title, x=0.01, y=0.92, ha="left", color=INK_2, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--datadir", required=True)
    ap.add_argument("--truthdir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    cfg = load_scenario(args.config)
    global PHRED_THRESHOLD
    PHRED_THRESHOLD = int(cfg.get("pipeline_phred_threshold", PHRED_THRESHOLD))
    truth, data, out = resolve(args.truthdir), resolve(args.datadir), resolve(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    haps = {r.id: str(r.seq) for r in SeqIO.parse(truth / "haplotypes.fasta", "fasta")}
    samples = pd.read_csv(truth / "samples.tsv", sep="\t")
    exp = pd.read_csv(truth / "expected_frequencies.tsv", sep="\t")
    frags = pd.read_csv(truth / "fragment_counts.tsv", sep="\t")
    lq_positions = {int(r["position"]) for r in cfg.get("low_quality_sites", []) or []}
    inserts = read_inserts(cfg)

    rows, lq_rows, mism, bases = [], [], 0, 0
    for s in samples.itertuples():
        counts, lq, m, b = count_reads(data / s.bam, haps, lq_positions)
        mism, bases = mism + m, bases + b
        for amp in inserts:
            total = sum(counts[(amp, h)] for h in haps)
            for h in haps:
                rows.append(dict(sample=s.sample, location=s.location, date=s.date, amplicon=amp, haplotype=h,
                                 reads=counts[(amp, h)], observed_frequency=counts[(amp, h)] / total if total else 0))
        for pos, (cov, low) in lq.items():
            lq_rows.append(dict(sample=s.sample, location=s.location, date=s.date, position=pos,
                                reads_covering=cov, reads_low_quality=low,
                                fraction_low_quality=low / cov if cov else 0))

    df = pd.DataFrame(rows).merge(exp, on=["location", "date", "haplotype"])
    df = df.merge(frags.rename(columns={"n_fragments": "simulated_fragments"}),
                  on=["sample", "location", "date", "amplicon", "haplotype"])
    df.to_csv(out / "observed_abundance.tsv", sep="\t", index=False)
    lq_df = pd.DataFrame(lq_rows)
    lq_df.to_csv(out / "low_quality_check.tsv", sep="\t", index=False)

    for d in (df, exp):
        d["date"] = pd.to_datetime(d["date"])
        d["clade"] = d["haplotype"].map(clade_of)
    plot(df, exp, cfg, "haplotype", out / "abundance_over_time.png",
         "Lines: expected frequency · points: observed in BAM (mean over amplicons, bars = range)")
    plot(df, exp, cfg, "clade", out / "clade_over_time.png",
         "Lines: expected frequency · points: observed in BAM (mean over amplicons, bars = range)")

    lost = (df.simulated_fragments - df.reads).sum()
    dev = (df.observed_frequency - df.expected_frequency).abs()
    print(f"samples: {len(samples)}  read pairs kept: {df.reads.sum()} / simulated {df.simulated_fragments.sum()} "
          f"(lost {lost})")
    print(f"observed vs expected frequency: mean |diff| {dev.mean():.3f}, max {dev.max():.3f}")
    print(f"read error rate vs true version: {mism / bases:.5f}")
    if len(lq_df):
        hit = lq_df[lq_df.reads_low_quality > 0]
        print(f"low-quality edge case: {len(hit)} samples with Q<={PHRED_THRESHOLD} bases at the site, "
              f"fraction {hit.fraction_low_quality.min():.2f}-{hit.fraction_low_quality.max():.2f}; "
              f"other samples: {int(lq_df[lq_df.reads_low_quality == 0].shape[0])} with none")
    print(f"plots: {out / 'abundance_over_time.png'}, {out / 'clade_over_time.png'}")


if __name__ == "__main__":
    main()
