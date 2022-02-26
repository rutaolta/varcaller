rule bwa_index:
    input:
        assembly=get_fasta
    output:
        amb=out_index_dir_path / "{sample}.amb",
        ann=out_index_dir_path / "{sample}.ann",
        bwt=out_index_dir_path / "{sample}.bwt",
        pac=out_index_dir_path / "{sample}.pac",
        sa=out_index_dir_path / "{sample}.sa"
    params:
    log:
        log=log_dir_path / "{sample}/bwa_index.log",
        cluster_log=cluster_log_dir_path / "{sample}/bwa_map.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/bwa_map.cluster.err"
    benchmark:
        benchmark=benchmark_dir_path / "{sample}/bwa_index.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bwa_index_threads"],
        time=config["bwa_index_time"],
        mem=config["bwa_index_mem_mb"]
    threads:
        config["bwa_index_threads"]
    shell:
        "bwa index {input.assembly} -p {out_index_dir_path}/{wildcards.sample} 2>{log.log}"
#TODO prefix


rule bwa_map:
    input:
        forward_reads=get_forward_fastq,
        reverse_reads=get_reverse_fastq,
        assembly=get_fasta,
        amb=rules.bwa_index.output.amb,
        ann=rules.bwa_index.output.ann,
        bwt=rules.bwa_index.output.bwt,
        pac=rules.bwa_index.output.pac,
        sa=rules.bwa_index.output.sa
    output:
        bam=out_alignment_dir_path / "{sample}/{reads}.sorted.mkdup.bam"
    params:
        bwa_threads=config["bwa_mem_threads"],
        fixmate_threads=config["samtools_fixmate_threads"],
        sort_threads=config["samtools_sort_threads"],
        markdup_threads=config["samtools_markdup_threads"],
        per_thread_sort_mem="%sG" % config["bwa_map_per_thread_mem_mb"],
        prefix=lambda wildcards, output: output["bam"][:-4]
    log:
        bwa_mem=log_dir_path / "{sample}/{reads}.bwa_mem.log",
        samtools_fixmate=log_dir_path / "{sample}/{reads}.samtools_fixmate.log",
        samtools_sort=log_dir_path / "{sample}/{reads}.samtools_sort.log",
        samtools_markdup=log_dir_path / "{sample}/{reads}.samtools_markdup.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.bwa_map.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.bwa_map.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/{reads}.bwa_map.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["bwa_mem_threads"] + config["samtools_fixmate_threads"] + config["samtools_sort_threads"] + config["samtools_markdup_threads"],
        time=config["bwa_map_time"],
        mem=config["bwa_mem_mem_mb"] + config["samtools_fixmate_mem_mb"] + config["bwa_map_per_thread_mem_mb"] * config["samtools_sort_threads"] * 1024 + config["samtools_markdup_mem_mb"]
    threads:
        config["bwa_mem_threads"] + config["samtools_fixmate_threads"] + config["samtools_sort_threads"] + config["samtools_markdup_threads"]
    shell:
        "bwa mem -t {params.bwa_threads} {out_index_dir_path}/{wildcards.sample} <(zcat -fc {input.forward_reads}) <(zcat -fc {input.reverse_reads}) "
        "-R  \'@RG\\tID:{wildcards.sample}\\tPU:x\\tSM:{wildcards.sample}\\tPL:Illumina\\tLB:x\' 2>{log.bwa_mem} | "
        "samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.samtools_fixmate} | "
        "samtools sort -T {params.prefix} -@ {params.sort_threads} -m {params.per_thread_sort_mem} 2>{log.samtools_sort} | "
        "samtools markdup -@ {params.markdup_threads} - {output.bam} 2>{log.samtools_markdup} "
#TODO platform or whole -R

rule index_bam:
    input:
        rules.bwa_map.output.bam
    output:
        bai=out_alignment_dir_path / "{sample}/{reads}.sorted.mkdup.bam.bai"
    log:
        std=log_dir_path / "{sample}/{reads}.index_bam.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.index_bam.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.index_bam.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/{reads}.index_bam.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["index_bam_threads"],
        time=config["index_bam_time"],
        mem=config["index_bam_mem_mb"],
    threads:
        config["index_bam_threads"]
    shell:
        "samtools index -@ {threads} {input} > {log.std} 2>&1"
