from snakemake.utils import min_version
from pathlib import Path
from os import walk

##### set minimum snakemake version #####
min_version("5.4.0")

##### setup config #####
configfile: "config/default.yaml"

# input dirs
assemblies_dir_path = Path(config["assemblies_dir"])
reads_dir_path = Path(config["reads_dir"])
# output dirs
out_alignment_dir_path = Path(config["out_alignment_dir"])
# technical dirs
scripts_dir_path = str(config["scripts_dir"])
# log dirs
log_dir_path = Path(config["log_dir"])
cluster_log_dir_path = Path(config["cluster_log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])

# def get_scaffolds(mypath):
#     _, _, filenames = next(walk(mypath))
#     return [splitext(filename)[0] for filename in filenames if filename.endswith('.fasta')]

# SAMPLES = get_scaffolds(reads_dir_path)

SAMPLES=config["reads"]

#### load rules #####
include: "workflow/rules/Alignment/alignment.smk"
include: "workflow/rules/Alignment/coverage.smk"

##### target rules #####
localrules: all

rule all:
    input:
        expand(out_alignment_dir_path / "{sample}/{sample}.coverage.per-base.bed.gz", sample=SAMPLES)
