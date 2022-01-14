rule bwa_index:
    input:
        assembly=assemblies_dir_path / (config["assembly"] + "." + config["assemblies_ext"] + ".gz")
    output:
        amb=temp(assemblies_dir_path / (config["assembly"] + "." + config["assemblies_ext"] + ".amb")),
        ann=temp(assemblies_dir_path / (config["assembly"] + "." + config["assemblies_ext"] + ".ann")),
        bwt=temp(assemblies_dir_path / (config["assembly"] + "." + config["assemblies_ext"] + ".bwt")),
        pac=temp(assemblies_dir_path / (config["assembly"] + "." + config["assemblies_ext"] + ".pac")),
        sa=temp(assemblies_dir_path / (config["assembly"] + "." + config["assemblies_ext"] + ".sa")),
        assembly=temp(assemblies_dir_path / (config["assembly"] + "." + config["assemblies_ext"]))
    params:
        prefix=assemblies_dir_path / (config["assembly"] + "." + config["assemblies_ext"])
    log:
        bwa_index_log=log_dir_path / (config["assembly"] + ".bwa_index.log"),
        cluster_log=cluster_log_dir_path / (config["assembly"] + ".bwa_map.cluster.log"),
        cluster_err=cluster_log_dir_path / (config["assembly"] + ".bwa_map.cluster.err")
    benchmark:
        benchmark_dir_path / config["assembly"] / "bwa_index.benchmark.txt"
    conda:
       "../../envs/conda.yaml"
    resources:
        cpus=config["bwa_index_threads"],
        time=config["bwa_index_time"],
        mem=config["bwa_index_mem_mb"]
    threads:
        config["bwa_index_threads"]
    shell:
        "gunzip -k {input.assembly}; "
        "bwa index {output.assembly} -p {params.prefix} 2>{log.bwa_index_log}"

rule bwa_map:
    input:
        forward_reads=reads_dir_path / ("{sample}_1." + config["reads_ext"] + ".gz"),
        reverse_reads=reads_dir_path / ("{sample}_2." + config["reads_ext"] + ".gz"),
        assembly=rules.bwa_index.output.assembly,
        amb=rules.bwa_index.output.amb,
        ann=rules.bwa_index.output.ann,
        bwt=rules.bwa_index.output.bwt,
        pac=rules.bwa_index.output.pac,
        sa=rules.bwa_index.output.sa
    output:
        bam=temp(out_alignment_dir_path / "{sample}/{sample}.sorted.mkdup.bam")
    params:
        bwa_threads=config["bwa_mem_threads"],
        fixmate_threads=config["samtools_fixmate_threads"],
        sort_threads=config["samtools_sort_threads"],
        markdup_threads=config["samtools_markdup_threads"],
        per_thread_sort_mem="%sG" % config["bwa_map_per_thread_mem_mb"],
        prefix=expand(out_alignment_dir_path / "{sample}/{sample}", sample=SAMPLES)
    log:
        bwa_mem=log_dir_path / "{sample}/bwa_mem.log",
        samtools_fixmate=log_dir_path / "{sample}/samtools_fixmate.log",
        samtools_sort=log_dir_path / "{sample}/samtools_sort.log",
        samtools_markdup=log_dir_path / "{sample}/samtools_markdup.log",
        cluster_log=cluster_log_dir_path / "{sample}.bwa_map.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.bwa_map.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/bwa_map.benchmark.txt"
    conda:
       "../../envs/conda.yaml"
    resources:
        cpus=config["bwa_mem_threads"] + config["samtools_fixmate_threads"] + config["samtools_sort_threads"] + config["samtools_markdup_threads"] + 1,
        time=config["bwa_map_time"],
        mem=config["bwa_mem_mem_mb"] + config["samtools_fixmate_mem_mb"] + config["bwa_map_per_thread_mem_mb"] * config["samtools_sort_threads"] * 1024 + config["samtools_markdup_mem_mb"]
    threads:
        config["bwa_mem_threads"] + config["samtools_fixmate_threads"] + config["samtools_sort_threads"] + config["samtools_markdup_threads"]
    shell:
        "bwa mem -t {params.bwa_threads} {input.assembly} <(zcat -c {input.forward_reads}) <(zcat -c {input.reverse_reads}) "
        "-R  \'@RG\\tID:{wildcards.sample}\\tPU:x\\tSM:{wildcards.sample}\\tPL:Illumina\\tLB:x\' 2>{log.bwa_mem} | "
        "samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.samtools_fixmate} | "
        "samtools sort -T {params.prefix} -@ {params.sort_threads} -m {params.per_thread_sort_mem} 2>{log.samtools_sort} | "
        "samtools markdup -@ {params.markdup_threads} - {output.bam} 2>{log.samtools_markdup} "
        "rm {input.assembly}"

rule index_bam:
    input:
        rules.bwa_map.output.bam
    output:
        bai=temp(out_alignment_dir_path / "{sample}/{sample}.sorted.mkdup.bam.bai")
    log:
        std=log_dir_path / "{sample}/index_bam.log",
        cluster_log=cluster_log_dir_path / "{sample}.index_bam.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.index_bam.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}/index_bam.benchmark.txt"
    conda:
        "../../envs/conda.yaml"
    resources:
        cpus=config["index_bam_threads"],
        time=config["index_bam_time"],
        mem=config["index_bam_mem_mb"],
    threads: 
        config["index_bam_threads"]
    shell:
        "samtools index -@ {threads} {input} > {log.std} 2>&1"
