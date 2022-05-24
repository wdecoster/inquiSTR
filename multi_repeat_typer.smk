#inquiSTR
import os
import glob
import pandas as pd


def process_targets(bed):
    targets = pd.read_table(bed, dtype='string')
    targets["locus"] = targets.locus.apply(
        lambda x: x.replace(' ', '-').replace('(', '').replace(')', ''))
    targets.set_index("locus", drop=False, inplace=True)
    targets['region'] = targets["chromosome"] + ':' + targets["start"] + '-' + targets["end"]
    return targets


def get_region(wildcards):
    return targets.loc[{wildcards.locus}, 'region'].values[0]


#IDS, = glob_wildcards("/home/wdecoster/p200/phased_bams/{id}.bam")
IDS = [os.path.basename(p).replace('.bam', '') for p in glob.glob("/home/wdecoster/p200/phased_bams/*.bam") if os.path.isfile(p)]
dir = os.path.dirname(workflow.basedir)
outdir = "/home/wdecoster/p200/workflow_results/multi_repeat_typer/"
bed = os.path.join(dir, "data/repeats_for_repeat_typing.bed")
ref = "/home/wdecoster/GRCh38_recommended/GRCh38.fa"
targets = process_targets(bed)


rule all:
    input:
        os.path.join(outdir, "all_repeats.tsv"),
        os.path.join(outdir, "plot_all_repeat.html"),

rule repeat_typer:
    input:
        "/home/wdecoster/p200/phased_bams/{id}.bam"
    output:
        os.path.join(outdir, "repeats/{id}_{locus}.fa")
    log:
        os.path.join(outdir, "logs/repeat_typer_{id}_{locus}.log")
    params:
        ref = ref,
        script = os.path.join(dir, "scripts/repeat_typer.py"),
        region = get_region
    conda:
        os.path.join(dir, "envs/repeat_typer.yml")
    shell:
        "python {params.script} \
        --bam {input} \
        --ref {params.ref} \
        --locus {wildcards.locus} \
        --region {params.region} > {output} 2> {log}"

rule combine_calls:
    input:
        expand(os.path.join(outdir, "repeats/{id}_{locus}.fa"),
               id=IDS, locus=targets.index)
    output:
        os.path.join(outdir, "all_repeats.tsv")
    log:
        os.path.join(outdir, "logs/combine_calls.log")
    params:
        script = os.path.join(dir, "scripts/combine_repeat_calls.py")
    conda:
        os.path.join(dir, "envs/combine_repeats.yml")
    shell:
        "python {params.script} --fas {input} --output {output} 2> {log}"

rule add_sample_info:
    input:
        os.path.join(outdir, "all_repeats.tsv")
    output:
        os.path.join(outdir, "all_repeats_sample_info.tsv")
    log:
        os.path.join(outdir, "logs/add_sample_info.log")
    params:
        script = os.path.join(dir, "scripts/add_sample_info_to_all_repeats.py")
    conda:
        os.path.join(dir, "envs/add_sample_info.yml")
    shell:
        "python {params.script} {input} > {output} 2> {log}"


rule plot_lengths:
    input:
        os.path.join(outdir, "all_repeats_sample_info.tsv")
    output:
        os.path.join(outdir, "plot_all_repeat.html")
    log:
        os.path.join(outdir, "logs/plot_lengths.log")
    params:
        script = os.path.join(dir, "scripts/plot_repeat_length_distribution.py")
    conda:
        os.path.join(dir, "envs/plot_repeats.yml")
    shell:
        "python {params.script} {input} > {output} 2> {log}"
