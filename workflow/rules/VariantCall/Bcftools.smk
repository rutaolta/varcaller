ruleorder: bcftools_filter_indel_snp > bcftools_vcf_subset > bcftools_filter > bcftools_varcall   


rule bcftools_varcall:
    input:
        assembly=FASTA,
        samples=expand(out_alignment_dir_path / "{sample}/{assembly}.{sample}.sorted.mkdup.bam", assembly=ASSEMBLY, sample=SAMPLES.sample_id),
        indexes=expand(out_alignment_dir_path / "{sample}/{assembly}.{sample}.sorted.mkdup.bam.bai", assembly=ASSEMBLY, sample=SAMPLES.sample_id)
    output:
        mpileup=varcall_dir_path / (ASSEMBLY + ".mpileup.vcf.gz"),
        call=varcall_dir_path / (ASSEMBLY + ".vcf.gz")
    params:
        adjustMQ=50,
        annotate_mpileup=config["bcftools_mpileup_annotate"],
        annotate_call=config["bcftools_call_annotate"],
        max_depth=config["bcftools_mpileup_max_depth"],
        min_MQ=config["bcftools_mpileup_min_MQ"],
        min_BQ=config["bcftools_mpileup_min_BQ"]
    log:
        mpileup=log_dir_path / (ASSEMBLY + ".bcftools_mpileup.log"),
        call=log_dir_path / (ASSEMBLY + ".bcftools_call.log"),
        cluster_log=cluster_log_dir_path / (ASSEMBLY + ".bcftools_varcall.cluster.log"),
        cluster_err=cluster_log_dir_path / (ASSEMBLY + ".bcftools_varcall.cluster.err")
    benchmark:
        benchmark_dir_path / (ASSEMBLY + ".bcftools_varcall.benchmark.txt")
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
        mpileup=rules.bcftools_varcall.output.mpileup,
        call=rules.bcftools_varcall.output.call
    output:
        filt_vcf=varcall_dir_path / (ASSEMBLY + ".filt.vcf.gz")
    params:
        soft_filter=config["bcftools_filter_soft_filter"],
        exclude=config["bcftools_filter_exclude"],
    log:
        std=log_dir_path / (ASSEMBLY + ".bcftools_filter.log"),
        cluster_log=cluster_log_dir_path / (ASSEMBLY + ".bcftools_filter.cluster.log"),
        cluster_err=cluster_log_dir_path / (ASSEMBLY + ".bcftools_filter.cluster.err")
    benchmark:
        benchmark_dir_path / (ASSEMBLY + ".bcftools_filter.benchmark.txt")
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_threads"],
        mem=config["bcftools_filter_mem_mb"],
        time=config["bcftools_filter_time"]
    threads:
        config["bcftools_filter_threads"]
    shell:
        "bcftools filter -Oz -s {params.soft_filter} --exclude '{params.exclude}' {input.call} > {output.filt_vcf} 2> {log.std}; "


checkpoint bcftools_vcf_subset:
    input:
        rules.bcftools_filter.output.filt_vcf
    output:
        dir=directory(vcf_subset_dir_path)
    params:
        variants=vcf_subset_dir_path / (ASSEMBLY + ".csv")
    log:
        std=log_dir_path / (ASSEMBLY + ".bcftools_vcf_subset.log"),
        cluster_log=cluster_log_dir_path / (ASSEMBLY + ".bcftools_vcf_subset.cluster.log"),
        cluster_err=cluster_log_dir_path / (ASSEMBLY + ".bcftools_vcf_subset.cluster.err")
    benchmark:
        benchmark_dir_path / (ASSEMBLY + ".bcftools_vcf_subset.benchmark.txt")
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_vcf_subset_threads"],
        mem=config["bcftools_vcf_subset_mem_mb"],
        time=config["bcftools_vcf_subset_time"]
    threads:
        config["bcftools_vcf_subset_threads"]
    shell:
        "for SUBSET in `bcftools query -l {input}`; "
        "do mkdir -p {output.dir}/$SUBSET; "
        "echo $SUBSET >> {params.variants}; "
        "bcftools view -Oz -s $SUBSET {input} > {output.dir}/$SUBSET/{ASSEMBLY}.vcf.gz 2> {log.std}; "
        "done; "


rule bcftools_filter_indel_snp:
    input:
        subvcf=vcf_subset_dir_path / "{subset}/{assembly}.vcf.gz"
    output:
        indel=vcf_subset_dir_path / "{subset}/{assembly}.indel.vcf.gz",
        snp=vcf_subset_dir_path / "{subset}/{assembly}.snp.vcf.gz"
    params:
        type_indel=config["bcftools_filter_indel_type"],
        type_snp=config["bcftools_filter_snp_type"]
    log:
        std=log_dir_path / "{subset}/{assembly}.bcftools_filter_indel_snp.log",
        cluster_log=cluster_log_dir_path / "{subset}/{assembly}.bcftools_filter_indel_snp.cluster.log",
        cluster_err=cluster_log_dir_path / "{subset}/{assembly}.bcftools_filter_indel_snp.cluster.err"
    benchmark:
        benchmark_dir_path / "{subset}/{assembly}.bcftools_filter_indel_snp.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_indel_snp_threads"],
        mem=config["bcftools_filter_indel_snp_mem_mb"],
        time=config["bcftools_filter_indel_snp_time"]
    threads:
        config["bcftools_filter_indel_snp_threads"]
    shell:
        "bcftools  filter -i {params.type_indel} -Oz {input.subvcf} > {output.indel} 2> {log.std}; "
        "bcftools  filter -i {params.type_snp} -Oz {input.subvcf} > {output.snp} 2> {log.std}; "


