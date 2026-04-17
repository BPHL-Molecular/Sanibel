process seroba {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/seroba", mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        path("seroba_output/pred.csv")
        val meta, emit: done

    script:
    """
    seroba runSerotyping /seroba/database/ ${reads[0]} ${reads[1]} seroba_output
    """
}
