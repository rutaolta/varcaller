rule mosdepth:
    input:
        bam=rules.bwa_map.output.bam,
        bai=rules.index_bam.output.bai
    output:
        out=alignment_dir_path / "{sample}/{assembly}.{sample}.coverage.per-base.bed.gz"
    params:
        min_mapping_quality=config["min_mapping_quality"],
        prefix=lambda w: alignment_dir_path / ("{sample}/{assembly}.{sample}.coverage".format(assembly = w.assembly, sample=w.sample))
    log:
        std=log_dir_path / "{sample}/{assembly}.{sample}.mosdepth.log",
        cluster_log=cluster_log_dir_path / "{sample}/{assembly}.{sample}.mosdepth.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{assembly}.{sample}.mosdepth.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{assembly}.{sample}.mosdepth.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["mosdepth_threads"],
        time=config["mosdepth_time"],
        mem=config["mosdepth_mem_mb"],
    threads: config["mosdepth_threads"]
    shell:
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {params.prefix} {input.bam} > {log.std} 2>&1"


rule coverage_whole_genome_stats:
    input:
        rules.mosdepth.output.out #alignment_dir_path / "{sample}/{assembly}.{sample}.coverage.per-base.bed.gz"
    output:
        alignment_dir_path / "{sample}/{assembly}.{sample}.coverage_whole_genome_stats.csv"
    params:
        prefix=lambda w: alignment_dir_path / ("{sample}/{assembly}.{sample}.coverage".format(assembly = w.assembly, sample=w.sample))
    log:
        std=log_dir_path / "{sample}/{assembly}.{sample}.coverage_whole_genome_stats.log",
        cluster_log=cluster_log_dir_path / "{sample}/{assembly}.{sample}.coverage_whole_genome_stats.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{assembly}.{sample}.coverage_whole_genome_stats.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{assembly}.{sample}.coverage_whole_genome_stats.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["coverage_stats_threads"],
        time=config["coverage_stats_time"],
        mem=config["coverage_stats_mem_mb"],
    threads:
        config["coverage_stats_threads"]
    shell:
        "coverage_statistics.py -i {input} --tool-name mosdepth -g -o {params.prefix} > {log.std} 2>&1"


rule coverage_window_stats:
    input:
        rules.mosdepth.output.out #alignment_dir_path / "{sample}/{assembly}.{sample}.coverage.per-base.bed.gz"
    output:
        alignment_dir_path / ("{sample}/{assembly}.{sample}.coverage_{size}_windows_stats.csv")
    params:
        window_size=lambda w: "{size}".format(size=w.size),
        prefix=lambda w: alignment_dir_path / ("{sample}/{assembly}.{sample}.coverage".format(assembly = w.assembly, sample=w.sample)),
    log:
        std=log_dir_path / "{sample}/{assembly}.{sample}.{size}.coverage_window_stats.log",
        cluster_log=cluster_log_dir_path / "{sample}/{assembly}.{sample}.{size}.coverage_window_stats.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{assembly}.{sample}.{size}.coverage_window_stats.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{assembly}.{sample}.{size}.coverage_window_stats.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["coverage_stats_threads"],
        time=config["coverage_stats_time"],
        mem=config["coverage_stats_mem_mb"],
    threads:
        config["coverage_stats_threads"]
    shell:
        "coverage_statistics.py -i {input} --tool-name mosdepth -n -f {params.window_size} -o {params.prefix} 2>&1"


rule coverage_visualization:
    input:
        whole_stats=rules.coverage_whole_genome_stats.output, #alignment_dir_path / "{sample}/{assembly}.{sample}.coverage_whole_genome_stats.csv",
        window_stats=rules.coverage_window_stats.output, #alignment_dir_path / ("{sample}/{assembly}.{sample}.{size}.coverage_windows_stats.csv"),
        length=assembly_stats_dir_path / "{assembly}.len",
        whitelist=assembly_stats_dir_path / "{assembly}.whitelist",
        syn=assembly_stats_dir_path / "{assembly}.syn",
        renamelist=assembly_stats_dir_path / "{assembly}.renamelist"
    output:
        png=alignment_dir_path / ("{sample}/{assembly}.{sample}.{size}.track.jet.png"),
        svg=alignment_dir_path / ("{sample}/{assembly}.{sample}.{size}.track.jet.svg")
    params:
        prefix=lambda w: alignment_dir_path / ("{sample}/{assembly}.{sample}.{size}.track").format(assembly=w.assembly, sample=w.sample, size=w.size),
        window_size=lambda w: "{size}".format(size=w.size),
        label=lambda w: "{sample}/{assembly}.{sample}.{size}".format(assembly=w.assembly, sample=w.sample, size=w.size),
    log:
        std=log_dir_path / "{sample}/{assembly}.{sample}.{size}.coverage_visualization.log",
        cluster_log=cluster_log_dir_path / "{sample}/{assembly}.{sample}.{size}.coverage_visualization.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{assembly}.{sample}.{size}.coverage_visualization.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{assembly}.{sample}.{size}.coverage_visualization.benchmark.txt"
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


