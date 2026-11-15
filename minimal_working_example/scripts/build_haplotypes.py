"""Build the haplotype ("version") sequences in reference coordinates.

Source FASTAs are aligned to the reference and projected onto it (substitutions only;
reference bases fill positions the source does not cover). Derived versions copy their
parent and add the listed substitutions.

Outputs (in --outdir):
  haplotypes.fasta            one full-length sequence per version, reference coordinates
  haplotype_mutations.tsv     every difference to the reference, per version
  pairwise_distances.tsv      SNP distances between versions within the amplicon inserts
"""
import argparse
import itertools
import sys

from Bio import Align, SeqIO

from common import load_scenario, read_inserts, read_reference, resolve

TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def project_onto_reference(ref, query, name):
    aligner = Align.PairwiseAligner(mode="global", match_score=2, mismatch_score=-3,
                                    open_gap_score=-6, extend_gap_score=-2, end_gap_score=0)
    aln = aligner.align(ref, query)[0]
    hap = list(ref)
    ref_blocks, query_blocks = aln.aligned
    if len(ref_blocks) > 1:
        print(f"WARNING {name}: alignment has {len(ref_blocks) - 1} gap(s); "
              "insertions are dropped, deleted positions keep the reference base", file=sys.stderr)
    for (rs, re_), (qs, _) in zip(ref_blocks, query_blocks):
        for i in range(re_ - rs):
            hap[rs + i] = query[qs + i]
    covered = (int(ref_blocks[0][0]), int(ref_blocks[-1][1]))
    return hap, covered


def parse_mutation(spec, hap):
    spec = str(spec)
    pos = int(spec.rstrip("ACGT"))
    alt = spec[len(str(pos)):] or TRANSITION[hap[pos - 1]]
    return pos, alt


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    cfg = load_scenario(args.config)
    ref_id, ref = read_reference(cfg)
    inserts = read_inserts(cfg)
    outdir = resolve(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def in_inserts(pos1):
        return [a for a, (s, e) in inserts.items() if s < pos1 <= e]

    haps = {}
    for name, h in cfg["haplotypes"].items():
        if "source" in h:
            query = str(next(SeqIO.parse(resolve(h["source"]), "fasta")).seq).upper()
            haps[name], covered = project_onto_reference(ref, query, name)
            for amp, (s, e) in inserts.items():
                if s < covered[0] or e > covered[1]:
                    sys.exit(f"{name}: source does not cover {amp} ({s}-{e}); covered {covered}")
        else:
            haps[name] = list(haps[h["parent"]])
        for spec in h.get("mutations", []):
            pos, alt = parse_mutation(spec, haps[name])
            amps = in_inserts(pos)
            if not amps:
                sys.exit(f"{name}: mutation at {pos} is outside all amplicon inserts")
            if haps[name][pos - 1] == alt:
                sys.exit(f"{name}: position {pos} already is {alt}")
            others = {haps[o][pos - 1] for o in haps if o != name}
            if len(others) > 1:
                print(f"WARNING {name}: position {pos} already varies between versions", file=sys.stderr)
            haps[name][pos - 1] = alt
            print(f"{name}: {ref[pos - 1]}{pos}{alt} (in {', '.join(amps)})")

    with open(outdir / "haplotypes.fasta", "w") as f:
        for name, seq in haps.items():
            f.write(f">{name}\n{''.join(seq)}\n")

    with open(outdir / "haplotype_mutations.tsv", "w") as f:
        f.write("haplotype\tposition\tref\talt\tamplicons\n")
        for name, seq in haps.items():
            for i, (r, a) in enumerate(zip(ref, seq)):
                if r != a:
                    f.write(f"{name}\t{i + 1}\t{r}\t{a}\t{','.join(in_inserts(i + 1))}\n")

    insert_pos = sorted({p for s, e in inserts.values() for p in range(s, e)})
    with open(outdir / "pairwise_distances.tsv", "w") as f:
        f.write("haplotype_a\thaplotype_b\tsnps_in_inserts\n")
        print("\nSNP distances within amplicon inserts:")
        for a, b in itertools.combinations(haps, 2):
            d = sum(haps[a][p] != haps[b][p] for p in insert_pos)
            f.write(f"{a}\t{b}\t{d}\n")
            print(f"  {a:12s} {b:12s} {d}")


if __name__ == "__main__":
    main()
