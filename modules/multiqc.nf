process multiqc {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/multiqc" }, mode: 'copy'

    input:
        tuple val(meta), path(reports)
    output:
        tuple val(meta), path("multiqc_report.html"), emit: report

    script:
    """
    multiqc .
    """
}

process multiqc_global {
    tag "multiqc_global"
    publishDir { "${params.output}" }, mode: 'copy'

    input:
        path summary

    output:
        path("interactive_report.html"), emit: report

    script:
    """
    multiqc ${params.output} \\
        --filename interactive_report.html \\
        --interactive \\
        --ignore "*/multiqc/*" \\
        --ignore "*interactive_report*"
    """
}
