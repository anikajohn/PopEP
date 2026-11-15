"""Shared helpers for the simulation scripts."""
import datetime as dt
import math
from pathlib import Path

import yaml
from Bio import SeqIO

ROOT = Path(__file__).resolve().parent.parent


def load_scenario(path):
    with open(path) as f:
        cfg = yaml.safe_load(f)
    return cfg


def resolve(path):
    """Resolve a scenario path relative to the minimal_example/ folder."""
    p = Path(path)
    return p if p.is_absolute() else ROOT / p


def read_reference(cfg):
    rec = next(SeqIO.parse(resolve(cfg["reference"]), "fasta"))
    return rec.id, str(rec.seq).upper()


def read_inserts(cfg):
    """Return {amplicon: (start0, end0)} with BED (0-based, end-exclusive) coordinates."""
    inserts = {}
    for amp, bed in cfg["amplicons"].items():
        with open(resolve(bed)) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    fields = line.rstrip("\n").split("\t")
                    inserts[amp] = (int(fields[1]), int(fields[2]))
                    break
    return inserts


def sampling_dates(cfg):
    s = cfg["sampling"]
    start, end = _date(s["start"]), _date(s["end"])
    dates, d = [], start
    while d <= end:
        dates.append(d)
        d += dt.timedelta(days=int(s["every_days"]))
    return dates


def _date(x):
    return x if isinstance(x, dt.date) else dt.date.fromisoformat(str(x))


def expected_frequencies(cfg, location, date):
    """Expected frequency of every haplotype at one location and date."""
    start = _date(cfg["locations"][location]["d31_start"])
    weights = {}
    for name, h in cfg["haplotypes"].items():
        if "growth_per_day" in h:
            days = (date - start).days
            weights[name] = 0.0 if days < 0 else h["seed_weight"] * math.exp(h["growth_per_day"] * days)
        else:
            weights[name] = float(h["weight"])
    total = sum(weights.values())
    return {k: v / total for k, v in weights.items()}


def sample_name(cfg, location, date):
    code = cfg["locations"][location]["code"]
    return f"{cfg['sampling']['sample_prefix']}_{code}_{date:%Y_%m_%d}"


def batch_name(date):
    return f"{date:%Y%m%d}_SIM"


def clade_of(name):
    return name.rsplit("_v", 1)[0]
