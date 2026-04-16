process summary_report {
    tag "summary"
    publishDir "${params.output}", mode: 'copy'

    input:
        path assembly_stats_files   // *_assembly_stats.txt  — also acts as BMGAP2+species gate
        path read_metrics_files     // *_readMetrics.txt
        path prokka_txt_files       // {sample_id}.txt  (prokka CDS count)
        path mlst_files             // {sample_id}.mlst
        path kraken_reports         // {sample_id}.report
        path pmga_files             // {sample_id}sta.txt
        path neisseria_txt
        path hinfluenzae_txt

    output:
        path "sum_report*.txt", emit: summary

    script:
    def outdir = params.output
    """
    summary_report.py \
        --outdir   "${outdir}" \
        --neisseria_txt   ${neisseria_txt} \
        --hinfluenzae_txt ${hinfluenzae_txt}
    """
}
