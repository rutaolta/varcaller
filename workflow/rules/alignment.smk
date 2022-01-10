rule bwa_index:
    input:
        assembly=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"] + ".gz")
    output:
        amb=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"] + ".amb"),
        ann=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"] + ".ann"),
        bwt=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"] + ".bwt"),
        pac=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"] + ".pac"),
        sa=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"] + ".sa")
    params:
        prefix=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"]),
        assembly=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"])
    log:
        bwa_index_log=log_dir_path / (config["sample"] + ".bwa_index.log"),
        cluster_log=cluster_log_dir_path / (config["sample"] + ".bwa_map.cluster.log"),
        cluster_err=cluster_log_dir_path / (config["sample"] + ".bwa_map.cluster.err")
    conda:
       "../envs/conda.yaml"
    resources:
        cpus=config["bwa_index_threads"],
        time=config["bwa_index_time"],
        mem=config["bwa_index_mem_mb"]
    threads:
        config["bwa_index_threads"]
    shell:
        "gunzip -k {input.assembly}; "
        "bwa index {params.assembly} -p {params.prefix} 2>{log.bwa_index_log}"

rule bwa_map:
    input:
        forward_reads=reads_dir_path / ("{sample}_1." + config["reads_ext"] + ".gz"),
        reverse_reads=reads_dir_path / ("{sample}_2." + config["reads_ext"] + ".gz"),
        assembly=assemblies_dir_path / (config["sample"] + "." + config["assemblies_ext"]),
        ann=rules.bwa_index.output.ann
    output:
        bam=temp(out_alignment_dir_path / "{sample}/{sample}.sorted.mkdup.bam")
    params:
        per_thread_sort_mem="%sG" % config["bwa_map_per_thread_mem_mb"],
        prefix=out_alignment_dir_path / "{sample_id}/{sample_id}"
    log:
        bwa_mem=log_dir_path / "{sample}/bwa_mem.log",
        samtools_fixmate=log_dir_path / "{sample}/samtools_fixmate.log",
        samtools_sort=log_dir_path / "{sample}/samtools_sort.log",
        samtools_markdup=log_dir_path / "{sample}/samtools_markdup.log",
        cluster_log=cluster_log_dir_path / "{sample}.bwa_map.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.bwa_map.cluster.err"
    conda:
       "../envs/conda.yaml"
    resources:
        cpus=config["bwa_map_threads"],
        time=config["bwa_map_time"],
        mem=config["bwa_map_mem_mb"]
    threads:
        config["bwa_map_threads"]
    shell:
        "bwa mem -t {threads} {input.assembly} <(zcat -c {input.forward_reads}) <(zcat -c {input.reverse_reads}) "
        "-R  \'@RG\\tID:{wildcards.sample}\\tPU:x\\tSM:{wildcards.sample}\\tPL:Illumina\\tLB:x\' 2>{log.bwa_mem} | "
        "samtools fixmate -@ {threads} -m - -  2>{log.samtools_fixmate} | "
        "samtools sort -T {params.prefix} -@ {threads} -m {params.per_thread_sort_mem} 2>{log.samtools_sort} | "
        "samtools markdup -@ {threads} - {output.bam} 2>{log.samtools_markdup} "
        "rm {input.assembly}"

rule index_bam:
    input:
        rules.bwa_map.output.bam
    output:
        out_alignment_dir_path / "{sample}/{sample}.sorted.mkdup.bam.bai"
    log:
        std=log_dir_path / "{sample}/index_bam.log",
        cluster_log=cluster_log_dir_path / "{sample}.index_bam.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.index_bam.cluster.err"
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
