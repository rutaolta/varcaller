import os

ruleorder: bcftools_varcall > bcftools_filter > bcftools_vcf_subset > bcftools_filter_indel_snp

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
        rules.bcftools_filter.output[0] #varcall_dir_path / "{sample}/{reads}.filt.vcf.gz"
    output:
        dir=directory(vcf_subset_dir_path / "{sample}/{reads}"),
        variants=vcf_subset_dir_path / "{sample}/{reads}.csv"
    params:
        prefix=lambda w, output: os.path.splitext(output.variants)[0]
    log:
        std=log_dir_path / "{sample}/{reads}.bcftools_vcf_subset.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.bcftools_vcf_subset.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.bcftools_vcf_subset.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/{reads}.bcftools_vcf_subset.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_vcf_subset_threads"],
        mem=config["bcftools_vcf_subset_mem_mb"],
        time=config["bcftools_vcf_subset_time"]
    threads:
        config["bcftools_vcf_subset_threads"]
    shell:
        "mkdir {output.dir}; "
        "for VARIANT in `bcftools query -l {input}`; "
        "do echo $VARIANT >> {output.variants}; "
        "bcftools view -Oz -s $VARIANT {input} > {output.dir}/$VARIANT.vcf.gz; "
        "done; "


rule bcftools_filter_indel_snp:
    input:
        subvcf=vcf_subset_dir_path / "{sample}/{reads}/{variant}.vcf.gz"
    output:
        indel=vcf_subset_dir_path / "{sample}/{reads}/{variant}.indel.vcf.gz",
        snp=vcf_subset_dir_path / "{sample}/{reads}/{variant}.snp.vcf.gz"
    params:
        type_indel=config["bcftools_filter_indel_type"],
        type_snp=config["bcftools_filter_snp_type"]
    log:
        std=log_dir_path / "{sample}/{reads}/{variant}.bcftools_filter_indel_snp.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}/{variant}.bcftools_filter_indel_snp.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}/{variant}.bcftools_filter_indel_snp.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/{reads}/{variant}.bcftools_filter_indel_snp.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_indel_snp_threads"],
        mem=config["bcftools_filter_indel_snp_mem_mb"],
        time=config["bcftools_filter_indel_snp_time"]
    threads:
        config["bcftools_filter_indel_snp_threads"]
    shell:
        "bcftools  filter -i {params.type_indel} -Oz {input.subvcf} > {output.indel}; "
        "bcftools  filter -i {params.type_snp} -Oz {input.subvcf} > {output.snp} "

