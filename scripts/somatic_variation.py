from os.path import basename
from Bio import SeqIO
from argparse import ArgumentParser
import pandas as pd
import numpy as np


def main():
    args = get_args()
    df = parse_repeats(args.fas)
    df.to_csv(args.output, sep="\t", na_rep="NA")


def parse_repeats(fastafiles):
    output = []
    for fas in fastafiles:
        name = basename(fas).replace(".fa", "")
        lengths = {
            "1": [],
            "2": [],
            "Unp": [],
        }
        locus = name.split("_")[-1]
        name = name.replace(f"_{locus}", "")
        for rec in SeqIO.parse(fas, "fasta"):
            features = {
                r.split("=")[0]: r.split("=")[1]
                for r in rec.description.split(" ")
                if not r == rec.id
            }
            flag = features["flag"] if "flag" in features.keys() else ""
            if features["type"] in ["contraction", "insert"] and flag != "aligned_softclip":
                lengths[features["phase"]].append(int(features["length"]))
        output.append(
            (
                name,
                locus,
                determine_variation(lengths["1"]),
                determine_variation(lengths["2"]),
            )
        )
    return pd.DataFrame(
        output,
        columns=[
            "individual",
            "repeat locus",
            "HP1",
            "HP2",
        ],
    ).set_index("individual")


def determine_variation(spanning_lengths):
    """if there are enough calls for this haplotype,
    return the median for spanning reads
    or the maximum of clipped reads"""
    if len(spanning_lengths) > 1:
        median = np.median(spanning_lengths)
        return np.mean([abs(median - i) for i in spanning_lengths])
    else:
        return np.nan


def get_args():
    parser = ArgumentParser(description="look at somatic variation in calling repeats")
    parser.add_argument("--fas", help="fasta files from repeat_typer.py", nargs="+", required=True)
    parser.add_argument("-o", "--output", help="output table", default="somatic_variation.tsv")
    return parser.parse_args()


if __name__ == "__main__":
    main()
