"""Compare PopEP raw_read_extraction read versions with the simulation truth.

Within each amplicon, versions with identical sequence are merged into one group
(e.g. "C.1.9.3_v1+C.1.9.3_v2" in amplicons that do not contain position 400). Each read version
(sequence + count) is assigned to the group with the fewest mismatches over its non-N positions;
ties are "ambiguous" (e.g. an N exactly at the only site separating two groups).

Usage (popep_sim env):
  python pipeline_test/compare_pipeline_output.py --config config/scenario.yaml \
      --truthdir results/truth --pipeline results/pipeline/raw_read_extraction --outdir results/pipeline_check
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from common import load_scenario, read_inserts, resolve  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--truthdir", required=True)
    ap.add_argument("--pipeline", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    cfg = load_scenario(args.config)
    truth, pipe, out = resolve(args.truthdir), resolve(args.pipeline), resolve(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    inserts = read_inserts(cfg)
    haps = {r.id: str(r.seq) for r in SeqIO.parse(truth / "haplotypes.fasta", "fasta")}
    lq_sites = [int(r["position"]) for r in cfg.get("low_quality_sites", []) or []]
    frags = pd.read_csv(truth / "fragment_counts.tsv", sep="\t")
    samples = pd.read_csv(truth / "samples.tsv", sep="\t")

    rows, group_of = [], {}
    for location in cfg["locations"]:
        for amp, (s, e) in inserts.items():
            fa = pipe / location / amp / "combined_amp_reads_versions.fasta"
            tsv = pipe / location / amp / "combined_amp_read_versions_counts.tsv"
            counts = pd.read_csv(tsv, sep="\t")
            counts = counts[counts.read_name != "read_name"]  # header repeated per concatenated file
            counts = dict(zip(counts.read_name, counts["count"].astype(int)))
            groups = {}
            for h, seq in haps.items():
                groups.setdefault(seq[s:e], []).append(h)
            amp_haps = {"+".join(hs): seq for seq, hs in groups.items()}
            group_of[amp] = {h: g for g in amp_haps for h in g.split("+")}
            for rec in SeqIO.parse(fa, "fasta"):
                seq = str(rec.seq)
                sample = rec.id.split("/")[0]
                dist = {h: sum(a != b for a, b in zip(seq, hs) if a != "N") for h, hs in amp_haps.items()}
                best = min(dist.values())
                hits = [h for h, d in dist.items() if d == best]
                n_at_site = any(s < p <= e and seq[p - 1 - s] == "N" for p in lq_sites)
                rows.append(dict(sample=sample, location=location, amplicon=amp, read_version=rec.id,
                                 count=counts[rec.id], assigned=hits[0] if len(hits) == 1 else "ambiguous",
                                 mismatches=best, n_total=seq.count("N"), n_at_lq_site=n_at_site))
    rv = pd.DataFrame(rows)
    rv.to_csv(out / "read_versions_assigned.tsv", sep="\t", index=False)

    obs = rv.groupby(["sample", "location", "amplicon", "assigned"])["count"].sum().rename("pipeline_reads")
    frags["assigned"] = [group_of[a][h] for a, h in zip(frags.amplicon, frags.haplotype)]
    truth_counts = (frags.groupby(["sample", "location", "amplicon", "assigned"]).n_fragments.sum()
                    .rename("simulated").reset_index())
    cmp = truth_counts.merge(obs.reset_index(), on=["sample", "location", "amplicon", "assigned"], how="outer")
    cmp[["simulated", "pipeline_reads"]] = cmp[["simulated", "pipeline_reads"]].fillna(0).astype(int)
    cmp.to_csv(out / "pipeline_vs_truth.tsv", sep="\t", index=False)

    # Summary
    n_samples = rv["sample"].nunique()
    print(f"samples in pipeline output: {n_samples} / {len(samples)}")
    tot = rv["count"].sum()
    print(f"reads in read versions: {tot} / simulated {frags.n_fragments.sum()}")
    print(f"reads in read versions with 0 mismatches to a true version: {rv.loc[rv.mismatches == 0, 'count'].sum() / tot:.4f}")
    amb = rv.loc[rv.assigned == "ambiguous", "count"].sum()
    print(f"reads ambiguous between versions: {amb}")
    known = cmp[cmp.assigned != "ambiguous"]
    diff = (known.pipeline_reads - known.simulated).abs()
    print(f"per sample x amplicon x version group |pipeline - simulated| reads: mean {diff.mean():.2f}, max {diff.max()}")
    print("version groups per amplicon:")
    for amp, g in group_of.items():
        print(f"  {amp}: {sorted(set(g.values()))}")

    # Low-quality edge case
    lq_reads = pd.read_csv(truth / "low_quality_reads.tsv", sep="\t")
    exp_n = lq_reads.groupby("sample").size().rename("expected_reads_with_N")
    got_n = rv[rv.n_at_lq_site].groupby("sample")["count"].sum().rename("pipeline_reads_with_N")
    n_cmp = pd.concat([exp_n, got_n], axis=1).fillna(0).astype(int)
    n_cmp.to_csv(out / "low_quality_N_check.tsv", sep="\t")
    print("\nN at low-quality site (expected vs pipeline):")
    print(n_cmp.to_string())
    other_n = rv[(~rv.n_at_lq_site) & (rv.n_total > 0)]["count"].sum()
    print(f"reads with N elsewhere: {other_n}")


if __name__ == "__main__":
    main()
