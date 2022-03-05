rule draw_variant_window_densities:
    input:
        hetero_indel=vcf_subset_dir_path / "{subset}/{assembly}.indel.hetero.vcf.gz",
        homo_indel=vcf_subset_dir_path / "{subset}/{assembly}.indel.homo.vcf.gz",
        hetero_snp=vcf_subset_dir_path / "{subset}/{assembly}.snp.hetero.vcf.gz",
        homo_snp=vcf_subset_dir_path / "{subset}/{assembly}.snp.homo.vcf.gz",
        whitelist=assembly_stats_dir_path / "{assembly}.whitelist",
        syn=assembly_stats_dir_path / "{assembly}.syn",
        renamelist=assembly_stats_dir_path / "{assembly}.renamelist"
    output:
        hetero_indel_png=vcf_subset_dir_path / "{subset}/{assembly}.indel.hetero.{size_and_step}.jet.png",
        hetero_indel_svg=vcf_subset_dir_path / "{subset}/{assembly}.indel.hetero.{size_and_step}.jet.svg",
        hetero_indel_scaffolds_absent_in_reference=vcf_subset_dir_path / "{subset}/{assembly}.indel.hetero.{size_and_step}.scaffolds_absent_in_reference.ids",
        hetero_indel_scaffolds_absent_in_vcf=vcf_subset_dir_path / "{subset}/{assembly}.indel.hetero.{size_and_step}.scaffolds_absent_in_vcf.ids",
        hetero_indel_short_scaffolds=vcf_subset_dir_path / "{subset}/{assembly}.indel.hetero.{size_and_step}.short_scaffolds.ids",
        hetero_indel_variant_counts_stats=vcf_subset_dir_path / "{subset}/{assembly}.indel.hetero.{size_and_step}.variant_counts.stats",
        hetero_indel_variant_counts_tsv=vcf_subset_dir_path / "{subset}/{assembly}.indel.hetero.{size_and_step}.variant_counts.tsv"
    params:
        hetero_indel_prefix=lambda w: vcf_subset_dir_path / ("{subset}/{assembly}.indel.hetero").format(assembly=w.assembly, subset=w.subset),
        homo_indel_prefix=lambda w: vcf_subset_dir_path / ("{subset}/{assembly}.indel.homo").format(assembly=w.assembly, subset=w.subset),
        hetero_snp_prefix=lambda w: vcf_subset_dir_path / ("{subset}/{assembly}.snp.hetero").format(assembly=w.assembly, subset=w.subset),
        homo_snp_prefix=lambda w: vcf_subset_dir_path / ("{subset}/{assembly}.snp.homo").format(assembly=w.assembly, subset=w.subset),
        density_thresholds=config["density_thresholds"],
        subplots_adjust_left=config["subplots_adjust_left"],
        size_and_step=lambda w: "{size_and_step}".format(size_and_step=w.size_and_step),
        syn_file_key_column=config["syn_file_key_column"],
        syn_file_value_column=config["syn_file_value_column"]
    log:
        std=log_dir_path / "{subset}/{assembly}.{size_and_step}.draw_variant_window_densities.log",
        cluster_log=cluster_log_dir_path / "{subset}/{assembly}.{size_and_step}.draw_variant_window_densities.cluster.log",
        cluster_err=cluster_log_dir_path / "{subset}/{assembly}.{size_and_step}.draw_variant_window_densities.cluster.err"
    benchmark:
        benchmark_dir_path / "{subset}/{assembly}.{size_and_step}.draw_variant_window_densities.benchmark.txt"
    conda:
        "../../../%s" % config["conda_config"]
    resources:
        cpus=config["draw_variant_window_densities_threads"],
        mem=config["draw_variant_window_densities_mem_mb"],
        time=config["draw_variant_window_densities_time"]
    threads:
        config["draw_variant_window_densities_threads"]
    shell:
        "draw_variant_window_densities.py -i {input.hetero_indel} -o SNPs.hetero.100000.vcf --density_thresholds '{params.density_thresholds}' --hide_track_label "
        "--subplots_adjust_left '{params.subplots_adjust_left}' -l 'Variant density' -w '{params.size_and_step}' -s '{params.size_and_step}' -a {input.whitelist} -z {input.renamelist} "
        " --scaffold_syn_file {input.syn} --syn_file_key_column '{params.syn_file_key_column}' --syn_file_value_column '{params.syn_file_value_column}' --colormap jet 2> {log.std}; "
