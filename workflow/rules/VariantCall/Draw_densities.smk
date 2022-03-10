rule draw_variant_window_densities:
    input:
        vcf=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.vcf.gz",
        whitelist=assembly_stats_dir_path / "{assembly}.whitelist",
        syn=assembly_stats_dir_path / "{assembly}.syn",
        renamelist=assembly_stats_dir_path / "{assembly}.renamelist"
    output:
        png=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.jet.png",
        svg=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.jet.svg",
        scaffolds_absent_in_reference=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.scaffolds_absent_in_reference.ids",
        scaffolds_absent_in_vcf=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.scaffolds_absent_in_vcf.ids",
        short_scaffolds=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.short_scaffolds.ids",
        variant_counts_stats=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.variant_counts.stats",
        variant_counts_tsv=vcf_subset_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.variant_counts.tsv"
    params:
        prefix=lambda w: vcf_subset_dir_path / ("{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}").format(assembly=w.assembly, subset=w.subset, ploidy=w.ploidy, var_type=w.var_type, zygosity=w.zygosity, size_and_step=w.size_and_step),
        density_thresholds=config["density_thresholds"],
        subplots_adjust_left=config["subplots_adjust_left"],
        syn_file_key_column=config["syn_file_key_column"],
        syn_file_value_column=config["syn_file_value_column"]
    log:
        std=log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.draw_variant_window_densities.log",
        cluster_log=cluster_log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.draw_variant_window_densities.cluster.log",
        cluster_err=cluster_log_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.draw_variant_window_densities.cluster.err"
    benchmark:
        benchmark_dir_path / "{subset}/{assembly}.{ploidy}.{var_type}.{zygosity}.{size_and_step}.draw_variant_window_densities.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["draw_variant_window_densities_threads"],
        mem=config["draw_variant_window_densities_mem_mb"],
        time=config["draw_variant_window_densities_time"]
    threads:
        config["draw_variant_window_densities_threads"]
    shell:
        "draw_variant_window_densities.py -i {input.vcf} -o {params.prefix} --density_thresholds '{params.density_thresholds}' --hide_track_label "
        "--subplots_adjust_left '{params.subplots_adjust_left}' -l 'Variant density' -w '{wildcards.size_and_step}' -s '{wildcards.size_and_step}' -a {input.whitelist} -z {input.renamelist} "
        " --scaffold_syn_file {input.syn} --syn_file_key_column '{params.syn_file_key_column}' --syn_file_value_column '{params.syn_file_value_column}' --colormap jet 2> {log.std}; "
