def gather_from_checkpoint(wildcards):
    checkpoint_output = checkpoints.bcftools_vcf_subset.get(**wildcards).output["dir"]
    return expand(vcf_subset_dir_path / "{sample}/{reads}/{variant}.indel.vcf.gz",
           sample=wildcards.sample,
           reads=wildcards.reads,
           variant=glob_wildcards(os.path.join(checkpoint_output, "{variant}.indel.vcf.gz")).variant)


rule fake_rule:
    input:
        gather_from_checkpoint
    output:
        Path("data_output/result") / "{sample}/{reads}.txt"
    log:
        std=log_dir_path / "{sample}/{reads}.bcftools_filter_indel_snp.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.bcftools_filter_indel_snp.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.bcftools_filter_indel_snp.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/{reads}.bcftools_filter_indel_snp.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_indel_snp_threads"],
        mem=config["bcftools_filter_indel_snp_mem_mb"],
        time=config["bcftools_filter_indel_snp_time"]
    threads:
        config["bcftools_filter_indel_snp_threads"]
    shell:
        "touch {output} "


