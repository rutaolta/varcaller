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
        "../../../envs/conda.yaml"
    resources:
        cpus=config["mosdepth_threads"],
        time=config["mosdepth_time"],
        mem=config["mosdepth_mem_mb"],
    threads: config["mosdepth_threads"]
    shell:
        "mosdepth -t {threads} --mapq {params.min_mapping_quality} {out_alignment_dir_path}/{wildcards.sample}/{wildcards.reads}.coverage {input.bam} > {log.std} 2>&1"
#TODO prefix