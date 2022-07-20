import fdb
import sys
from argparse import ArgumentParser
import pandas as pd


class Gentli(object):
    def __init__(self, user, pw):
        self.usr = user
        self.pw = pw
        self.role = "RR_SC"
        self.con = fdb.connect(
            dsn="molgenvz.cde.ua.ac.be:/home/firebird/gentli.fdb",
            user=self.usr,
            password=self.pw,
            role=self.role,
        )
        self.cur = self.con.cursor()

    def run_to_individual(self, run_id):
        reply = self.cur.execute(
            'select "ont_individual" ' 'from "RR_SC:full_ont" ' 'where "ngs_seq_id" = ?', (run_id,)
        ).fetchall()
        return reply[0][0]


def main():
    args = get_args()
    cohort = parse_cohort()
    gentli_login = {
        line.split("=")[0]: line.rstrip().split("=")[1]
        for line in open("/home/wdecoster/p200/scripts/secrets.txt")
    }
    g = Gentli(gentli_login["user"], gentli_login["pw"])
    repeats = open(args.all_repeats)
    header = next(repeats).rstrip()
    print(f"{header}\tgentli_id\tphenotype")
    for line in repeats:
        individual = g.run_to_individual(line.split("\t")[0].split("_")[-1].split("u")[0])
        try:
            path = cohort.loc[individual, "path"]
        except KeyError:
            path = "Unavailable"
        print(f"{line.rstrip()}\t{individual}\t{path}")


def parse_cohort():
    df = pd.read_excel("/home/wdecoster/cohorts/Individuals.xlsx")
    df["path"] = df["PathDx1"].str.cat(df["TDP43_Type"].fillna(""))
    controls = ["Normal", "AGD", "SC", "CVD", "UP", "CTE", "PA", "PART"]
    others = [
        i
        for i in df["path"].unique()
        if not i in ["aFTLD-U", "FTLD-TDPA", "FTLD-TDPB", "FTLD-TDPC"] + controls
    ]
    return (
        df[["Gentli_ID", "path"]]
        .set_index("Gentli_ID")
        .replace(others, "other")
        .replace(controls, "control")
    )


def get_args():
    parser = ArgumentParser()
    parser.add_argument("all_repeats", help="file created by combining repeats of repeat_typer")
    return parser.parse_args()


if __name__ == "__main__":
    main()
