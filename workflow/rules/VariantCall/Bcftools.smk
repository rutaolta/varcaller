# ruleorder: bcftools_vcf_subset > create_subset_out_dirs > bcftools_filter_indel_snp > bcftools_varcall > bcftools_filter
ruleorder: bcftools_vcf_subset > bcftools_filter_hetero_homo > bcftools_filter_indel_snp > bcftools_varcall > bcftools_filter


rule bcftools_varcall:
    input:
        assembly=FASTA,
        samples=expand(alignment_dir_path / "{sample}/{assembly}.{sample}.sorted.mkdup.bam", assembly=ASSEMBLY, sample=SAMPLES.sample_id),
        indexes=expand(alignment_dir_path / "{sample}/{assembly}.{sample}.sorted.mkdup.bam.bai", assembly=ASSEMBLY, sample=SAMPLES.sample_id)
    output:
        mpileup=varcall_dir_path / (ASSEMBLY + "." + PLOIDY + ".mpileup.vcf.gz"),
        call=varcall_dir_path / (ASSEMBLY + "." + PLOIDY + ".vcf.gz")
    params:
        adjustMQ=50,
        annotate_mpileup=config["bcftools_mpileup_annotate"],
        annotate_call=config["bcftools_call_annotate"],
        max_depth=config["bcftools_mpileup_max_depth"],
        min_MQ=config["bcftools_mpileup_min_MQ"],
        min_BQ=config["bcftools_mpileup_min_BQ"],
        samples_file=assembly_stats_dir_path / (ASSEMBLY + ".samples.file") if ploidy_of_ChrX else '',
        regions_file=assembly_stats_dir_path / (ASSEMBLY + ".ploidy.file") if ploidy_of_ChrX else ''
    log:
        mpileup=log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_mpileup.log"),
        call=log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_call.log"),
        cluster_log=cluster_log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_varcall.cluster.log"),
        cluster_err=cluster_log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_varcall.cluster.err")
    benchmark:
        benchmark_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_varcall.benchmark.txt")
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
        "--adjust-MQ {params.adjustMQ} --annotate {params.annotate_mpileup} -Oz "
        "--samples-file {params.samples_file} --regions-file {params.regions_file} "
        "-f {input.assembly} {input.samples} 2> {log.mpileup} | "
        "tee {output.mpileup} | bcftools call -Oz -mv --annotate {params.annotate_call} > {output.call} 2> {log.call}"


rule bcftools_filter:
    input:
        mpileup=rules.bcftools_varcall.output.mpileup,
        call=rules.bcftools_varcall.output.call
    output:
        filt_vcf=varcall_dir_path / (ASSEMBLY + "." + PLOIDY + ".filt.vcf.gz")
    params:
        soft_filter=config["bcftools_filter_soft_filter"],
        exclude=config["bcftools_filter_exclude"],
    log:
        std=log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_filter.log"),
        cluster_log=cluster_log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_filter.cluster.log"),
        cluster_err=cluster_log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_filter.cluster.err")
    benchmark:
        benchmark_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_filter.benchmark.txt")
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
        variants=vcf_subset_dir_path / (ASSEMBLY + "." + PLOIDY + ".csv")
    log:
        std=log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_vcf_subset.log"),
        cluster_log=cluster_log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_vcf_subset.cluster.log"),
        cluster_err=cluster_log_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_vcf_subset.cluster.err")
    benchmark:
        benchmark_dir_path / (ASSEMBLY + "." + PLOIDY + ".bcftools_vcf_subset.benchmark.txt")
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
        "bcftools view -Oz -s $SUBSET {input} > {output.dir}/$SUBSET/{ASSEMBLY}.{PLOIDY}.vcf.gz 2> {log.std}; "
        "done; "


rule bcftools_filter_indel_snp:
    input:
        subvcf=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.vcf.gz"
    output:
        var_vcf=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.vcf.gz"
    params:
        variance_type=lambda w: config["bcftools_filter_variance_type"].format(var_type=w.var_type)
    log:
        std=log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.bcftools_filter_indel_snp.log",
        cluster_log=cluster_log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.bcftools_filter_indel_snp.cluster.log",
        cluster_err=cluster_log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.bcftools_filter_indel_snp.cluster.err"
    benchmark:
        benchmark_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.bcftools_filter_indel_snp.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_indel_snp_threads"],
        mem=config["bcftools_filter_indel_snp_mem_mb"],
        time=config["bcftools_filter_indel_snp_time"]
    threads:
        config["bcftools_filter_indel_snp_threads"]
    shell:
        "bcftools  filter -i '{params.variance_type}' -Oz {input.subvcf} > {output.var_vcf} 2> {log.std}; "


rule bcftools_filter_hetero_homo:
    input:
        var_vcf=rules.bcftools_filter_indel_snp.output.var_vcf
    output:
        var_zyg_vcf=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.vcf.gz"
    params:
        zygosity_type=lambda w: config["bcftools_filter_zygosity_type"].format(zygosity=w.zygosity)
    log:
        std=log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.bcftools_filter_hetero_homo.log",
        cluster_log=cluster_log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.bcftools_filter_hetero_homo.cluster.log",
        cluster_err=cluster_log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.bcftools_filter_hetero_homo.cluster.err"
    benchmark:
        benchmark_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.bcftools_filter_hetero_homo.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bcftools_filter_hetero_homo_threads"],
        mem=config["bcftools_filter_hetero_homo_mem_mb"],
        time=config["bcftools_filter_hetero_homo_time"]
    threads:
        config["bcftools_filter_hetero_homo_threads"]
    shell:
        "bcftools filter -i '{params.zygosity_type}' -Oz {input.var_vcf} > {output.var_zyg_vcf} 2> {log.std}; "
