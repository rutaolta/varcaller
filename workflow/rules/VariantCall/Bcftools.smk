ruleorder: bcftools_varcall > bcftools_filter

rule bcftools_varcall:
    input:
        assembly=get_fasta,
        samples=rules.bwa_map.output.bam,
        indexes=rules.index_bam.output.bai
    output:
        mpileup=varcall_dir_path / "{sample}/{reads}.mpileup.vcf.gz",
        call=varcall_dir_path / "{sample}/{reads}.vcf.gz"
    params:
        adjustMQ=50,
        annotate_mpileup=config["bcftools_mpileup_annotate"],
        annotate_call=config["bcftools_call_annotate"],
        max_depth=config["bcftools_mpileup_max_depth"],
        min_MQ=config["bcftools_mpileup_min_MQ"],
        min_BQ=config["bcftools_mpileup_min_BQ"]
    log:
        mpileup=log_dir_path / "{sample}/{reads}.bcftools_mpileup.log",
        call=log_dir_path / "{sample}/{reads}.bcftools_call.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.bcftools_varcall.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.bcftools_varcall.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/{reads}.bcftools_varcall.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_varcall_threads"],
        mem=config["bcftools_varcall_mem_mb"],
        time=config["bcftools_varcall_time"]
    threads:
        config["bcftools_varcall_threads"]
    shell:
        "bcftools mpileup --threads {threads} -d {params.max_depth} -q {params.min_MQ} -Q {params.min_BQ} "
        "--adjust-MQ {params.adjustMQ} --annotate {params.annotate_mpileup} -Oz -f {input.assembly} {input.samples} 2> {log.mpileup} | "
        "tee {output.mpileup} | bcftools call -Oz -mv --annotate {params.annotate_call} > {output.call} 2> {log.call}"
#TODO error: Call: unrecognized option '--annotate'
#TODO error: Mpileup: Could not parse tag "SCR" in "AD,INFO/AD,ADF,INFO/ADF,ADR,INFO/ADR,DP,SP,SCR,INFO/SCR"


rule bcftools_filter:
    input:
        rules.bcftools_varcall.output.call
    output:
        varcall_dir_path / "{sample}/{reads}.filt.vcf.gz"
    params:
        soft_filter=config["bcftools_filter_soft_filter"],
        exclude=config["bcftools_filter_exclude"],
    log:
        std=log_dir_path / "{sample}/{reads}.bcftools_filter.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.bcftools_filter.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.bcftools_filter.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/{reads}.bcftools_filter.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_threads"],
        mem=config["bcftools_filter_mem_mb"],
        time=config["bcftools_filter_time"]
    threads:
        config["bcftools_filter_threads"]
    shell:
        "bcftools filter -Oz -s {params.soft_filter} --exclude '{params.exclude}' {input} > {output} 2> {log.std}; "

checkpoint bcftools_vcf_subset:
    input:
        varcall_dir_path / (config["samples"] + ".filt.vcf.gz")
    output:
        directory(varcall_dir_path / "vcf_subset")
    log:
        std=log_dir_path / "bcftools_vcf_subset.log",
        cluster_log=cluster_log_dir_path / "bcftools_vcf_subset.cluster.log",
        cluster_err=cluster_log_dir_path / "bcftools_vcf_subset.cluster.err"
    benchmark:
        benchmark_dir_path / "bcftools_vcf_subset.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_vcf_subset_threads"],
        mem=config["bcftools_vcf_subset_mem_mb"],
        time=config["bcftools_vcf_subset_time"]
    threads:
        config["bcftools_vcf_subset_threads"]
    shell:
        "for SAMPLE in `bcftools query -l {input}`; "
        "do mkdir -p {output}/${{SAMPLE}}; "
        "bcftools view -Oz -s ${{SAMPLE}} {input} > {output}/${{SAMPLE}}/${{SAMPLE}}.vcf.gz; "
        "done"

ruleorder: bcftools_vcf_subset > create_out_dirs > bcftools_filter_indel_snp

rule bcftools_filter_indel_snp:
    input:
        subvcf=expand(vcf_subset_dir_path / "{SAMPLE}/{SAMPLE}.vcf.gz", SAMPLE=config["samples"])
    output:
        directory(vcf_subset_dir_path)
    params:
        type_indel=config["bcftools_filter_indel_type"],
        type_snp=config["bcftools_filter_snp_type"]
    log:
        std=log_dir_path / "bcftools_filter_indel_snp.log",
        cluster_log=cluster_log_dir_path / "bcftools_filter_indel_snp.cluster.log",
        cluster_err=cluster_log_dir_path / "bcftools_filter_indel_snp.cluster.err"
    benchmark:
        benchmark_dir_path / "bcftools_filter_indel_snp.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_indel_snp_threads"],
        mem=config["bcftools_filter_indel_snp_mem_mb"],
        time=config["bcftools_filter_indel_snp_time"]
    threads:
        config["bcftools_filter_indel_snp_threads"]
    shell:
        "bcftools  filter -i {params.type_indel} -Oz {input.subvcf} > {output}/{{SAMPLE}}.indel.vcf.gz; "
        "&& "
        "bcftools  filter -i {params.type_snp} -Oz {input.subvcf} > {output}/{{SAMPLE}}.snp.vcf.gz"
