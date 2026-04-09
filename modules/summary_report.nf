process summary_report {
    tag "summary"
    publishDir "${params.output}", mode: 'copy'

    input:
        path row_files

    output:
        path "sum_report*.txt", emit: summary

    script:
    """
    summary_report.py
    """
}
