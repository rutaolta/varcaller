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
varcall_dir_path = Path(config["varcall_dir"])
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

if "reads" not in config: # for parse samples IDs from directory
    config["reads"] = [d.name for d in reads_dir_path.iterdir() if d.is_dir()]
SAMPLES=config["reads"]

#### load rules #####
include: "workflow/rules/Alignment/alignment.smk"
include: "workflow/rules/Alignment/coverage.smk"
include: "workflow/rules/VariantCall/Bcftools.smk"

##### target rules #####
localrules: all

rule all:
    input:
        # Mosdepth:
        expand(out_alignment_dir_path / "{sample}/{sample}.coverage.per-base.bed.gz", sample=SAMPLES),
        # Variant call:
        expand(varcall_dir_path / (config["assembly"] + ".mpileup.vcf.gz")),
        expand(varcall_dir_path / (config["assembly"] + ".vcf.gz"))
        expand(varcall_dir_path / (config["assembly"] + ".filt.vcf.gz"))
