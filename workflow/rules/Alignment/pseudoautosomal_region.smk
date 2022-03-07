rule pseudoautosomal_region:
    input:
        whole_stats=rules.coverage_whole_genome_stats.output, #out_alignment_dir_path / "{sample}/{assembly}.{sample}.coverage_whole_genome_stats.csv",
        window_stats=out_alignment_dir_path / ("{sample}/{assembly}.{sample}.coverage_{par_size}_windows_stats.csv")
    output:
        bed=out_alignment_dir_path / "{sample}/PAR/{assembly}.{sample}.{par_size}_pseudoreg.bed",
        chrscaf=temp(out_alignment_dir_path / "{sample}/PAR/{assembly}.{sample}.{par_size}_chrscaf.csv")
    params:
        outdir=out_alignment_dir_path / "{sample}/PAR",
        prefix=lambda w: out_alignment_dir_path / ("{sample}/PAR/{assembly}.{sample}.{par_size}").format(assembly=w.assembly, sample=w.sample, par_size=PAR_SIZE),
        scaffold_name=config["sex_scaffold_name"],
        par_window_size=PAR_SIZE
    log:
        std=log_dir_path / "{sample}/{assembly}.{sample}.{par_size}.pseudoautosomal_region.log",
        cluster_log=cluster_log_dir_path / "{sample}/{assembly}.{sample}.{par_size}.pseudoautosomal_region.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}/{assembly}.{sample}.{par_size}.pseudoautosomal_region.cluster.err"
    benchmark:
         benchmark_dir_path / "{sample}/{assembly}.{sample}.{par_size}.pseudoautosomal_region.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["pseudoautosomal_region_threads"],
        time=config["pseudoautosomal_region_time"],
        mem=config["pseudoautosomal_region_mem_mb"]
    threads:
        config["pseudoautosomal_region_threads"]
    shell:
        "mkdir -p {params.outdir}; "
        "cat {input.window_stats} | awk '{{ if ($1 == \"'{params.scaffold_name}'\") print $0}}' > {output.chrscaf}; "
        "pseudoautosomal_region.py -f {params.par_window_size} -i {output.chrscaf} "
        "-s {params.scaffold_name} -o {params.prefix} "
        "-m $(cat {input.whole_stats} | sed -n 2p | awk '{{print $2}}') > {log.std} 2>&1 "


rule ploidy_file:
    input:
        beds=expand(out_alignment_dir_path / "{sample}/PAR/{assembly}.{sample}.{par_size}_pseudoreg.bed", assembly=ASSEMBLY, sample=SAMPLES.sample_id, par_size=PAR_SIZE),
        lenfile=rules.assembly_stats.output.lenfile
    output:
        varcall_dir_path / "ploidy.file"
    log:
        std=log_dir_path / "ploidy_file.log",
        cluster_log=cluster_log_dir_path / "ploidy_file.cluster.log",
        cluster_err=cluster_log_dir_path / "ploidy_file.cluster.err"
    benchmark:
         benchmark_dir_path / "ploidy_file.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["pseudoautosomal_region_threads"],
        time=config["pseudoautosomal_region_time"],
        mem=config["pseudoautosomal_region_mem_mb"]
    threads:
        config["pseudoautosomal_region_threads"]
    shell:
        "workflow/scripts/ploidy_file_generator.py --input {input.beds} --lenfile {input.lenfile} -o {output} > {log.std} 2>&1"