eval "$(/home/wdecoster/miniconda3/bin/conda shell.bash hook)"
conda activate snakemake


out=/home/wdecoster/p200/workflow_results/multi_repeat_typer/all_repeats.tsv
if [ -f $out ]; then
 rm $out
fi

si=/home/wdecoster/p200/workflow_results/multi_repeat_typer/all_repeats_sample_info.tsv
if [ -f $si ]; then
 rm $si
fi


plot=/home/wdecoster/p200/workflow_results/multi_repeat_typer/plot_all_repeat.html
if [ -f $plot ]; then
 rm $plot
fi

plot=/home/wdecoster/p200/workflow_results/multi_repeat_typer/somatic_variation.html``
if [ -f $plot ]; then
 rm $plot
fi


snakemake -s /home/wdecoster/p200/workflows/multi_repeat_typer.smk --cores 24 --use-conda
