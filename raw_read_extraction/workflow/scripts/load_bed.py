import pandas as pd

def load_bed(bedfile, sort=True):
    """
    Get coordinates of each amplicon region (start and stop coordinates).
    """
    amp_bed_raw = pd.read_table(
        bedfile,
        index_col=False,
        names=["ref", "start", "stop", "amp_num", "pool", "strand"],
        dtype={
            "ref": "str",
            "start": "uint32",
            "stop": "uint32",
            "amp_num": "str",
        },
    )
    amp_bed = amp_bed_raw.sort_values(by=["ref", "start"]) if sort else amp_bed_raw
    return amp_bed
