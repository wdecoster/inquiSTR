from argparse import ArgumentParser
import pandas as pd
from collections import Counter
from itertools import chain


def main():
    args = get_args()
    df = pd.read_table(args.outlier)
    c = Counter(
        list(
            chain.from_iterable(
                [
                    i.split(",")
                    for i in df["outliers"].str.replace("_H1", "").str.replace("_H2", "").values
                ]
            )
        )
    )
    for name, count in c.most_common():
        print(f"{name}\t{count}")


def get_args():
    parser = ArgumentParser("Quantify how often a sample is an outlier from inquiSTR outlier")
    parser.add_argument("outlier", help="file from inquiSTR outlier")
    return parser.parse_args()


if __name__ == "__main__":
    main()
