import os
import glob

IDS = [os.path.basename(p).replace('.bam', '') for p in glob.glob("/home/AD/wdecoster/p200/phased_bams/*.bam") if os.path.isfile(p)]
outdir = "/home/AD/wdecoster/p200"

rule all:
    input:
        expand(os.path.join(outdir, "inquistr-call/{id}.inq.gz"), id=IDS),
        os.path.join(outdir, "inquistr-combined/combined.inq.gz"),
        os.path.join(outdir, "inquistr-combined/outlier.inq.gz"),
        


rule inquistr_call:
    input:
        "/home/AD/wdecoster/p200/phased_bams/{id}.bam"
    output:
        os.path.join(outdir, "inquistr-call/{id}.inq.gz")
    log:
        os.path.join(outdir, "logs/inquistr-call/{id}.log")
    params:
        bed = "/home/AD/wdecoster/p200/simple_repeat_merged.bed"
    threads: 1
    shell:
        "/home/AD/wdecoster/p200/inquiSTR call -t {threads} -R {params.bed} {input} | gzip > {output} 2> {log}"

rule inquistr_combine:
    input:
        expand(os.path.join(outdir, "inquistr-call/{id}.inq.gz"), id=IDS),
    output:
        os.path.join(outdir, "inquistr-combined/combined.inq.gz")
    log:
        os.path.join(outdir, "logs/inquistr-combine.log")
    shell:
        "/home/AD/wdecoster/p200/inquiSTR combine {input} | gzip > {output} 2> {log}"

rule inquistr_outlier:
    input:
        os.path.join(outdir, "inquistr-combined/combined.inq.gz")
    output:
        os.path.join(outdir, "inquistr-combined/outlier.inq.gz")
    log:
        os.path.join(outdir, "logs/inquistr-outlier.log")
    shell:
        "/home/AD/wdecoster/p200/inquiSTR outlier {input} | gzip > {output} 2> {log}"