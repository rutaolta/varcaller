from snakemake.utils import min_version
from pathlib import Path
from os import walk

##### set minimum snakemake version #####
min_version("5.4.0")

##### setup config #####
configfile: "config/config.yaml"

# output dirs
out_index_dir_path = Path(config["out_index_dir"])
out_alignment_dir_path = Path(config["out_alignment_dir"])
varcall_dir_path = Path(config["varcall_dir"])
vcf_subset_dir_path = Path(config["vcf_subset_dir"])
# technical dirs
scripts_dir_path = str(config["scripts_dir"])
# log dirs
log_dir_path = Path(config["log_dir"])
cluster_log_dir_path = Path(config["cluster_log_dir"])
benchmark_dir_path = Path(config["benchmark_dir"])

# for bcftools_vcf_subset
def expand_template_from_bcftools_vcf_subset(wildcards, template):
    checkpoint_output = checkpoints.bcftools_vcf_subset.get(**wildcards).output[0]
    SAMPLE = glob_wildcards(os.path.join(checkpoint_output, "{SAMPLE}/{SAMPLE}.vcf.gz")).SAMPLE
    return expand(str(template), SAMPLE=SAMPLE)

#### load rules #####
include: "workflow/rules/common.smk"

include: "workflow/rules/Alignment/alignment.smk"
include: "workflow/rules/Alignment/coverage.smk"
include: "workflow/rules/VariantCall/Bcftools.smk"

##### target rules #####
localrules: all, create_out_dirs
ruleorder: create_out_dirs > bcftools_filter_indel_snp

rule all:
    input:
        # Mosdepth:
        expand(out_alignment_dir_path / "{sample}/{reads}.coverage.per-base.bed.gz", zip, sample=SAMPLES.assembly_name, reads=SAMPLES.forward_read_name),
        # Variant call:
        expand(varcall_dir_path / "{sample}/{reads}.filt.vcf.gz", zip, sample=SAMPLES.assembly_name, reads=SAMPLES.forward_read_name),
        lambda w: expand_template_from_bcftools_vcf_subset(w, vcf_subset_dir_path/ "{SAMPLE}.indel.vcf.gz"),
        lambda w: expand_template_from_bcftools_vcf_subset(w, vcf_subset_dir_path/ "{SAMPLE}.snp.vcf.gz")

rule create_out_dirs:
    output:
        vcf_subset_dirs=directory(expand(vcf_subset_dir_path / "{SAMPLE}", SAMPLE=SAMPLES))
    shell:
        "mkdir -p {output.vcf_subset_dirs}; "
