from snakemake.utils import min_version
from pathlib import Path
from os import walk

##### set minimum snakemake version #####
min_version("5.4.0")

##### setup config #####
configfile: "config/config.yaml"

# output dirs
assembly_stats_dir_path = Path(config["assembly_stats_dir"])
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


#### load rules #####
include: "workflow/rules/common.smk"

include: "workflow/rules/Alignment/alignment.smk"
include: "workflow/rules/Alignment/coverage.smk"
include: "workflow/rules/VariantCall/Bcftools.smk"
include: "workflow/rules/Preprocessing/assembly_stats.smk"


##### target rules #####
localrules: all, create_sample_cluster_log_dirs
ruleorder: create_sample_cluster_log_dirs > bcftools_vcf_subset

rule all:
    input:
        # alignment
        expand(out_alignment_dir_path / "{sample}/{assembly}.{sample}.sorted.mkdup.bam.bai", assembly=ASSEMBLY, sample=SAMPLES.sample_id),
        
        # variant calling
        lambda wildcards: aggregate_file_names(wildcards, str("{subset}/"+ASSEMBLY+".indel.vcf.gz")),
        lambda wildcards: aggregate_file_names(wildcards, str("{subset}/"+ASSEMBLY+".snp.vcf.gz")),

        # coverage visualization
        expand(out_alignment_dir_path / "{sample}/{assembly}.{sample}.{size}.track.jet.png", assembly=ASSEMBLY, sample=SAMPLES.sample_id, size=SIZE)


rule create_sample_cluster_log_dirs:
    output:
        directory(expand(cluster_log_dir_path / "{sample}", sample=SAMPLES.sample_id))
    shell:
        "mkdir -p {output}"


