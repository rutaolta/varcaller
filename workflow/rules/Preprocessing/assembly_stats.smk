rule samtools_faidx:
    input:
        assembly=FASTA
    output:
        assembly_stats_dir_path / "{assembly}.fai"
    log:
        std=log_dir_path / "{assembly}/samtools_faidx.log",
        cluster_log=cluster_log_dir_path / "{assembly}/samtools_faidx.cluster.log",
        cluster_err=cluster_log_dir_path / "{assembly}/samtools_faidx.cluster.err"
    benchmark:
        benchmark_dir_path / "{assembly}/samtools_faidx.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["samtools_faidx_threads"],
        mem=config["samtools_faidx_mem_mb"],
        time=config["samtools_faidx_time"]
    threads:
        config["samtools_faidx_threads"]
    shell:
        "samtools faidx {input.assembly} ;"
        "mv {input.assembly}.fai {output} "


rule assembly_stats:
    input:
        assembly_stats_dir_path / "{assembly}.fai"
    output:
        length=assembly_stats_dir_path / "{assembly}.len",
        whitelist=assembly_stats_dir_path / "{assembly}.whitelist",
        syn=assembly_stats_dir_path / "{assembly}.syn",
        renamelist=assembly_stats_dir_path / "{assembly}.renamelist"
    params:
        prefix=lambda wildcards, output: output["length"][:-4],
        chr_number=config["number_of_chromosomes"]
    log:
        std=log_dir_path / "{assembly}/assembly_stats.log",
        cluster_log=cluster_log_dir_path / "{assembly}/assembly_stats.cluster.log",
        cluster_err=cluster_log_dir_path / "{assembly}/assembly_stats.cluster.err"
    benchmark:
        benchmark_dir_path / "{assembly}/assembly_stats.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["samtools_faidx_threads"],
        mem=config["samtools_faidx_mem_mb"],
        time=config["samtools_faidx_time"]
    threads:
        config["samtools_faidx_threads"]
    shell:
        "cat {input} | awk '{{print $1\"\t\"$2}}' | sort -nr -k2 > {params.prefix}.len; "
        "cat {params.prefix}.len | awk '{{print $1}}' | head -n {params.chr_number} > {params.prefix}.whitelist; "
        "cat {params.prefix}.whitelist | awk '{{print $1\"\t\"$1}}' | head -n {params.chr_number} > {params.prefix}.syn; "
        "cat {params.prefix}.syn | awk '{{print $2}}' | head -n {params.chr_number} > {params.prefix}.renamelist "


