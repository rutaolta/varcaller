rule mosdepth:
    input:
        bam=rules.bwa_map.output.bam,
        bai=rules.index_bam.output.bai
    output:
        out=out_alignment_dir_path / "{sample}/{reads}.coverage.per-base.bed.gz"
    params:
        min_mapping_quality=config["min_mapping_quality"]
    log:
        std=log_dir_path / "{sample}/{reads}.mosdepth.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.mosdepth.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.mosdepth.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{reads}.mosdepth.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["mosdepth_threads"],
        time=config["mosdepth_time"],
        mem=config["mosdepth_mem_mb"],
    threads: config["mosdepth_threads"]
    shell:
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {out_alignment_dir_path}/{wildcards.sample}/{wildcards.reads}.coverage {input.bam} > {log.std} 2>&1"
#TODO prefix


rule coverage_whole_genome_stats:
    input:
        out_alignment_dir_path / "{sample}/{reads}.coverage.per-base.bed.gz"
    output:
        out_alignment_dir_path / "{sample}/{reads}_whole_genome_stats.csv"
    params:
        prefix=lambda wildcards, output: output[0][:-23]
    log:
        std=log_dir_path / "{sample}/{reads}.coverage_whole_genome_stats.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.coverage_whole_genome_stats.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.coverage_whole_genome_stats.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{reads}.coverage_whole_genome_stats.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["covarage_stats_threads"],
        time=config["covarage_stats_time"],
        mem=config["covarage_stats_mem_mb"],
    threads:
        config["covarage_stats_threads"]
    shell:
        "coverage_statistics.py -i {input} --tool-name mosdepth -g -o {params.prefix} > {log.std} 2>&1"


rule coverage_window_stats:
    input:
        out_alignment_dir_path / "{sample}/{reads}.coverage.per-base.bed.gz"
    output:
        out_alignment_dir_path / ("{sample}/{reads}_{size}_windows_stats.csv")
    params:
        window_size=lambda wildcards: "{size}".format(size=wildcards.size),
        prefix=lambda wildcards: out_alignment_dir_path / ("{sample}/{reads}".format(sample = wildcards.sample, reads=wildcards.reads)),
    log:
        std=log_dir_path / "{sample}/{reads}.{size}.coverage_window_stats.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.{size}.coverage_window_stats.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.{size}.coverage_window_stats.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{reads}.{size}.coverage_window_stats.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["covarage_stats_threads"],
        time=config["covarage_stats_time"],
        mem=config["covarage_stats_mem_mb"],
    threads:
        config["covarage_stats_threads"]
    shell:
        "coverage_statistics.py -i {input} --tool-name mosdepth -n -f {params.window_size} -o {params.prefix} 2>&1"


rule coverage_visualization:
    input:
        whole_stats=out_alignment_dir_path / "{sample}/{reads}_whole_genome_stats.csv",
        window_stats=out_alignment_dir_path / ("{sample}/{reads}_{size}_windows_stats.csv"),
        length=assembly_stats_dir_path / "{sample}.len",
        whitelist=assembly_stats_dir_path / "{sample}.whitelist",
        syn=assembly_stats_dir_path / "{sample}.syn",
        renamelist=assembly_stats_dir_path / "{sample}.renamelist"
    output:
        png=out_alignment_dir_path / ("{sample}/{reads}.{size}.track.jet.png"),
        svg=out_alignment_dir_path / ("{sample}/{reads}.{size}.track.jet.svg")
    params:
        prefix=lambda wildcards: out_alignment_dir_path / ("{sample}/{reads}.{size}.track").format(sample=wildcards.sample, reads=wildcards.reads, size=wildcards.size),
        window_size=lambda wildcards: "{size}".format(size=wildcards.size),
        label=lambda wildcards: "{sample}/{reads}.{size}".format(sample=wildcards.sample, reads=wildcards.reads, size=wildcards.size),
    log:
        std=log_dir_path / "{sample}/{reads}.{size}.coverage_visualization.log",
        cluster_log=cluster_log_dir_path / "{sample}/{reads}.{size}.coverage_visualization.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{reads}.{size}.coverage_visualization.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{reads}.{size}.coverage_visualization.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["coverage_visualization_threads"],
        time=config["coverage_visualization_time"],
        mem=config["coverage_visualization_mem_mb"]
    threads:
        config["coverage_visualization_threads"]
    shell:
        "draw_coverage.py --scaffold_column_name '#scaffold' --window_column_name 'frame' "
        "--coverage_column_name 'median' -i {input.window_stats} -o {params.prefix} "
        "--subplots_adjust_left 0.35 -l 'Coverage of {params.label}' "
        "-m $(cat {input.whole_stats} | sed -n 2p | awk '{{print $2}}') "
        "-w {params.window_size} -n {input.length} -a {input.whitelist} -z {input.renamelist} "
        "--scaffold_syn_file {input.syn} --syn_file_key_column 0 "
        "--syn_file_value_column 1 --colormap jet > {log.std} 2>&1; "