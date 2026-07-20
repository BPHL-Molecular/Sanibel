process summary_report {
    tag "summary"
    publishDir { "${params.output}" }, mode: 'copy'

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
        path skani_files
        path blast16s_files
        path amrfinder_files

    output:
        path "sum_report.txt",    emit: summary
        path "nm_sum_report.txt", emit: nm_summary, optional: true
        path "hi_sum_report.txt", emit: hi_summary, optional: true
        path "*_mqc.tsv",         emit: mqc_tables, optional: true

    script:
    def outdir = params.output
    """
    summary_report.py \
        --outdir   "${outdir}" \
        --neisseria_txt   ${neisseria_txt} \
        --hinfluenzae_txt ${hinfluenzae_txt}
    """
}
