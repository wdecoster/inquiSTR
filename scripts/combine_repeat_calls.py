from os.path import basename
from Bio import SeqIO
from argparse import ArgumentParser
import pandas as pd
import numpy as np


def main():
    args = get_args()
    repeats = (
        parse_repeats(args.fas, args.minsupport)
        .pivot_table(index="individual", columns="repeat locus", dropna=False)
        .sort_index(axis=1, level=1)
    )
    repeats.columns = [f"length_{locus}_{phase}" for phase, locus in repeats.columns]
    repeats.to_csv(args.output, sep="\t", na_rep="NA")


def parse_repeats(fastafiles, minsupport):
    output = []
    for fas in fastafiles:
        name = basename(fas).replace(".fa", "")
        lengths = {
            "1": {"spanning": [], "clipped": []},
            "2": {"spanning": [], "clipped": []},
            "Unp": {"spanning": [], "clipped": []},
        }
        locus = name.split("_")[-1]
        name = name.replace(f"_{locus}", "")
        for rec in SeqIO.parse(fas, "fasta"):
            features = {
                r.split("=")[0]: r.split("=")[1]
                for r in rec.description.split(" ")
                if not r == rec.id
            }
            evidence = (
                "spanning"
                if features["type"] in ["aligned_softclip", "contraction", "insert"]
                else "clipped"
            )
            lengths[features["phase"]][evidence].append(int(features["length"]))
        output.append(
            (
                name,
                locus,
                determine_genotype(lengths["1"], minsupport=minsupport),
                determine_genotype(lengths["2"], minsupport=minsupport),
                determine_genotype(lengths["Unp"], minsupport=minsupport),
            )
        )
    return pd.DataFrame(
        output, columns=["individual", "repeat locus", "HP1", "HP2", "Unp"]
    ).set_index("individual")


def determine_genotype(lengths_dict, minsupport=3):
    """if there are enough calls for this haplotype,
    return the median for spanning reads
    or the maximum of clipped reads"""
    if len(lengths_dict["spanning"] + lengths_dict["clipped"]) < minsupport:
        return np.nan
    elif lengths_dict["spanning"]:
        return np.median(lengths_dict["spanning"])
    else:
        return np.amax(lengths_dict["clipped"])


def get_args():
    parser = ArgumentParser(description="sum insertions")
    parser.add_argument(
        "--fas", help="fasta files from repeat_typer.py", nargs="+", required=True
    )
    parser.add_argument("-o", "--output", help="output table", default="summary.tsv")
    parser.add_argument(
        "-s",
        "--minsupport",
        help="minimal number of reads supporting a call",
        type=int,
        default=3,
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
