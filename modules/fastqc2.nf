process fastqc2 {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/fastqc2" }, mode: 'copy'

    input:
        tuple val(meta), path(clean_reads)
    output:
        tuple val(meta), path("*_clean_fastqc.{html,zip}"), emit: report

    script:
    """
    fastqc --threads ${task.cpus} ${clean_reads[0]} ${clean_reads[1]}
    """
}
