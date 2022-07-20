import pandas as pd
import plotly.express as px
from argparse import ArgumentParser


def main():
    args = get_args()
    df = pd.read_table(args.somatic_variation).melt(
        value_name="variation", var_name="haplotype", id_vars=["individual", "repeat locus"]
    )
    plot(df)


def plot(df):
    fig = px.strip(
        df,
        x="repeat locus",
        y="variation",
        hover_data=df.columns,
        labels={"Variation": "Mean absolute distance from median haplotype length [nucleotides]"},
    )
    print(fig.to_html(full_html=True, include_plotlyjs="cdn"))


def get_args():
    parser = ArgumentParser(description="make density plot of lengths of repeats")
    parser.add_argument("somatic_variation", help="TSV file somatic_variation.py")
    return parser.parse_args()


if __name__ == "__main__":
    main()
