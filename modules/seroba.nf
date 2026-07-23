process seroba {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/seroba" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        path("${meta.id}/pred.csv")
        val meta, emit: done

    script:
    """
    ln -s ${reads[0]} ${meta.id}_1.fastq.gz
    ln -s ${reads[1]} ${meta.id}_2.fastq.gz
    seroba runSerotyping /seroba/database/ ${meta.id}_1.fastq.gz ${meta.id}_2.fastq.gz ${meta.id}
    """
}
