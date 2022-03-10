from snakemake.utils import min_version
from pathlib import Path
from os import walk

##### set minimum snakemake version #####
min_version("5.4.0")

##### setup config #####
configfile: "config/config.yaml"

# output dirs
assembly_stats_dir_path = Path(config["assembly_stats_dir"])
index_dir_path = Path(config["out_index_dir"])
alignment_dir_path = Path(config["out_alignment_dir"])
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
include: "workflow/rules/Alignment/pseudoautosomal_region.smk"
include: "workflow/rules/VariantCall/Bcftools.smk"
include: "workflow/rules/VariantCall/Draw_densities.smk"
include: "workflow/rules/Preprocessing/assembly_stats.smk"


##### target rules #####
localrules: all, create_sample_cluster_log_dirs
ruleorder: create_sample_cluster_log_dirs > bcftools_vcf_subset

rule all:
    input:
        # create dirs for samples to avoid slurm error
        expand(cluster_log_dir_path / "{sample}", sample=SAMPLES.sample_id),

        # alignment
        expand(alignment_dir_path / "{sample}/{assembly}.{sample}.sorted.mkdup.bam.bai", assembly=ASSEMBLY, sample=SAMPLES.sample_id),

        # variant calling
        lambda wildcards: aggregate_file_names(str(vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.vcf.gz"), assembly=ASSEMBLY, ploidy=PLOIDY, var_type=VAR_TYPE),
        lambda wildcards: aggregate_file_names(str(vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.vcf.gz"), assembly=ASSEMBLY, ploidy=PLOIDY, var_type=VAR_TYPE, zygosity=ZYGOSITY),

        # draw densities
        lambda wildcards: aggregate_file_names(str(vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.jet.png"), assembly=ASSEMBLY, ploidy=PLOIDY, var_type=VAR_TYPE, zygosity=ZYGOSITY, size_and_step=SIZE_AND_STEP),

        # coverage visualization
        expand(alignment_dir_path / "{sample}/{assembly}.{sample}.{size}.track.jet.png", assembly=ASSEMBLY, sample=SAMPLES.sample_id, size=SIZE),

        # ploidy.file and samples.file
        expand(assembly_stats_dir_path / "{assembly}.ploidy.file", assembly=ASSEMBLY),
        expand(assembly_stats_dir_path / "{assembly}.samples.file", assembly=ASSEMBLY)

