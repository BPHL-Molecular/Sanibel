process multiqc {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/multiqc", mode: 'copy'

    input:
        tuple val(meta), path(reports)
    output:
        tuple val(meta), path("multiqc_report.html"), emit: report

    script:
    """
    multiqc .
    """
}
