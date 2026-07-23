process fastqc {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/fastqc" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("*_original_fastqc.{html,zip}"), emit: report

    script:
    def prefix = meta.id
    """
    fastqc --threads ${task.cpus} ${reads[0]} ${reads[1]}

    r1base=\$(basename ${reads[0]} .fastq.gz)
    r2base=\$(basename ${reads[1]} .fastq.gz)

    mv \${r1base}_fastqc.html ${prefix}_1_original_fastqc.html
    mv \${r1base}_fastqc.zip  ${prefix}_1_original_fastqc.zip
    mv \${r2base}_fastqc.html ${prefix}_2_original_fastqc.html
    mv \${r2base}_fastqc.zip  ${prefix}_2_original_fastqc.zip
    """
}
