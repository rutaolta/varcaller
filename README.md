# About

Pipeline gets assembly in fasta format, and reads in fq format. As an output it generates ... .

The ... you can find on the resulting plot in `data_output` folder.

# Configure Pipeline

`git clone https://github.com/rutaolta/reseqpipe.git`

`cd <pipeline_working_dir>`

It is recommended to create a fresh conda environment using mamba or conda.

```
mamba env create --name reseqpipe --file ./environment.yaml
# or:
# conda env create --name reseqpipe --file ./environment.yaml
```

Activate conda environment with snakemake:

`conda activate reseqpipe`

# Run

`snakemake --cores 32 --configfile config/default.yaml --forceall --use-conda --profile profile/slurm/ --printshellcmds --latency-wait 60`

