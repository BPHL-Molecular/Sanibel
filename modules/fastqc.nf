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
    mkdir -p raw
    ln -s ../${reads[0]} raw/${prefix}_R1.fastq.gz
    ln -s ../${reads[1]} raw/${prefix}_R2.fastq.gz

    fastqc --threads ${task.cpus} --outdir . raw/${prefix}_R1.fastq.gz raw/${prefix}_R2.fastq.gz

    mv ${prefix}_R1_fastqc.html ${prefix}_R1_original_fastqc.html
    mv ${prefix}_R1_fastqc.zip  ${prefix}_R1_original_fastqc.zip
    mv ${prefix}_R2_fastqc.html ${prefix}_R2_original_fastqc.html
    mv ${prefix}_R2_fastqc.zip  ${prefix}_R2_original_fastqc.zip
    """
}
