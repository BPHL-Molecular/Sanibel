process summary_report {
    tag "summary"
    publishDir "${params.output}", mode: 'copy'

    input:
        val  barrier
        path assembly_stats_files
        path read_metrics_files
        path prokka_txt_files
        path mlst_files
        path kraken_reports
        path pmga_files
        path neisseria_txt
        path hinfluenzae_txt
        path aggregate_files
        path skani_files
        path blast16s_files

    output:
        path "sum_report.txt",    emit: summary
        path "nm_sum_report.txt", emit: nm_summary, optional: true
        path "hi_sum_report.txt", emit: hi_summary, optional: true

    script:
    def outdir = params.output
    """
    summary_report.py \
        --outdir   "${outdir}" \
        --neisseria_txt   ${neisseria_txt} \
        --hinfluenzae_txt ${hinfluenzae_txt} \
        --min_ani      ${params.skani_routing_min_ani} \
        --min_af       ${params.skani_routing_min_af} \
        --min_coverage ${params.qc_min_coverage} \
        --warn_contigs ${params.qc_warn_contigs} \
        --fail_contigs ${params.qc_fail_contigs} \
        --min_n50      ${params.qc_min_n50}
    """
}
