from snakemake.utils import min_version
from pathlib import Path
from os import walk

##### set minimum snakemake version #####
min_version("5.4.0")

##### setup config #####
configfile: "config/config.yaml"

# input dirs
assemblies_dir_path = Path(config["assemblies_dir"])
# output dirs
out_alignment_dir_path = Path(config["out_alignment_dir"])
varcall_dir_path = Path(config["varcall_dir"])
# technical dirs
scripts_dir_path = str(config["scripts_dir"])
# log dirs
log_dir_path = Path(config["log_dir"])
cluster_log_dir_path = Path(config["cluster_log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])

#### load rules #####
include: "workflow/rules/common.smk"

include: "workflow/rules/Alignment/alignment.smk"
include: "workflow/rules/Alignment/coverage.smk"
include: "workflow/rules/VariantCall/Bcftools.smk"

##### target rules #####
localrules: all

rule all:
    input:
        # Mosdepth:
        expand(out_alignment_dir_path / "{sample}/{reads}.coverage.per-base.bed.gz", zip, sample=SAMPLES.assembly_name, reads=SAMPLES.forward_read_name),
        # Variant call:
        expand(varcall_dir_path / "{sample}/{reads}.filt.vcf.gz", zip, sample=SAMPLES.assembly_name, reads=SAMPLES.forward_read_name)
