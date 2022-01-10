# About

Pipeline gets assembly in fasta format, and reads in fq format as default (see `Required settings` section below). 
As an output it generates ... (pipeline is in progress).

All general results you can find in `data_output` folder.

# Configure Pipeline

`git clone https://github.com/rutaolta/reseqpipe.git`

`cd <pipeline_working_dir>`

It is recommended to create a fresh conda environment using `mamba` or `conda`.

```
mamba env create --name reseqpipe --file ./environment.yaml
# or:
# conda env create --name reseqpipe --file ./environment.yaml
```

Activate conda environment with snakemake:

`conda activate reseqpipe`

# Run

`snakemake --cores 32 --configfile config/default.yaml --forceall --use-conda --profile profile/slurm/ --printshellcmds --latency-wait 60`

# Required settings

Extensions of assembly and reads can be changed in `config/default.yaml` in "Extensions" section:

    - `assemblies_ext: "your desired extension"`
    
    - `reads_ext: "your desired extension"`

The name of file with reference assembly should be defined in `config/default.yaml` in "Parameters" section:

    - `sample: "your desired reference"`

The mapping quality for mosdepth can be changed in `config/default.yaml` in "Parameters" section:

    - `min_mapping_quality: your_desired_quality_integer` (20 as default)
