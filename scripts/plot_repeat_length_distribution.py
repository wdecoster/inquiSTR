import pandas as pd
import plotly.express as px
from argparse import ArgumentParser


def main():
    args = get_args()
    df = parse_repeats(args.repeats)
    plot(df)


def plot(df):
    fig = px.strip(
        df,
        x="locus",
        y="length",
        hover_data=df.columns,
        color="phenotype",
        labels={"length": "Expansion length [nucleotides]"},
    )
    print(fig.to_html(full_html=True, include_plotlyjs="cdn"))


def parse_repeats(repeatf):
    df = (
        pd.read_csv(repeatf, sep="\t", index_col="individual")
        .fillna(0)
        .drop_duplicates(subset="gentli_id")
    )
    df.drop(columns=[c for c in df.columns if c.endswith("_Unp")], inplace=True)
    lengths = {}
    for locus in columns_to_loci(df.columns):
        lengths[locus] = df[f"length_{locus}_HP1"].to_list() + df[f"length_{locus}_HP2"].to_list()
    return (
        pd.DataFrame(lengths, index=df.index.to_list() * 2)
        .melt(value_name="length", var_name="locus", ignore_index=False)
        .reset_index()
        .rename(columns={"index": "individual"})
        .sort_values(by="locus")
        .set_index("individual", drop=False)
        .join(df["phenotype"])
    )


def columns_to_loci(columns):
    return set(
        [
            c.replace("_HP1", "").replace("_HP2", "").replace("length_", "")
            for c in columns
            if c not in ["gentli_id", "phenotype"]
        ]
    )


def get_args():
    parser = ArgumentParser(description="make density plot of lengths of repeats")
    parser.add_argument("repeats", help="TSV file from combine_repeat_calls.py")
    return parser.parse_args()


if __name__ == "__main__":
    main()
