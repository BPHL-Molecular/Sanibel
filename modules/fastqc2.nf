process fastqc2 {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/fastqc2", mode: 'copy'

    input:
        tuple val(meta), path(clean_reads)
    output:
        tuple val(meta), path("*_clean_fastqc.{html,zip}"), emit: report

    script:
    def prefix = meta.id
    """
    fastqc --threads ${task.cpus} ${clean_reads[0]} ${clean_reads[1]}

    mv ${prefix}_1_fastqc.html ${prefix}_1_clean_fastqc.html
    mv ${prefix}_1_fastqc.zip  ${prefix}_1_clean_fastqc.zip
    mv ${prefix}_2_fastqc.html ${prefix}_2_clean_fastqc.html
    mv ${prefix}_2_fastqc.zip  ${prefix}_2_clean_fastqc.zip
    """
}
