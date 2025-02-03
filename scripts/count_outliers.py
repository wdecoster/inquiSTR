from argparse import ArgumentParser
import pandas as pd
from collections import Counter
from itertools import chain
import plotly.express as px
import numpy as np

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
    names = []
    counts = []
    for name, count in c.most_common():
        print(f"{name}\t{count}")
        names.append(name)
        counts.append(count)
    df = pd.DataFrame({"individual": names, "count": counts})
    if args.groups:
        groups = pd.read_table(args.groups, usecols=["individual", "group"], index_col="individual")
        df = df.join(groups, on="individual")

    fig = px.violin(
        df,
        x="group" if args.groups else None,
        color="group" if args.groups else None,
        y="count",
        title="Outlier loci<br>per individual",
        labels={"group": "", "count": "Number of outlier loci per individual"},
        box=True,
        points="all",
        hover_data=["individual"],
    )
    fig.update_traces(marker=dict(size=3), spanmode="hard")
    fig.update_layout(
        plot_bgcolor="white",
        font=dict(size=20),
        width=400,
        height=800,
        showlegend=False,
    )
    fig.update_xaxes(
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
    )
    fig.update_yaxes(
        showline=True,
        linewidth=2,
        linecolor="black",
        mirror=True,
        rangemode="tozero",
    )
    fig.write_html(args.output)


def get_args():
    parser = ArgumentParser("Quantify how often a sample is an outlier from inquiSTR outlier")
    parser.add_argument("outlier", help="file from inquiSTR outlier")
    parser.add_argument("--groups", help="file with group information")
    parser.add_argument("-o", "--output", help="output file", default="outliers.html")
    return parser.parse_args()


if __name__ == "__main__":
    main()
