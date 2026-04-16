process seroba {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/seroba", mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        path("seroba_output/pred.csv"), optional: true
        val meta, emit: done

    script:
    """
    if [[ "${meta.mash_species}" == "Streptococcus_pneumoniae" ]]; then
        seroba runSerotyping /seroba/database/ ${reads[0]} ${reads[1]} seroba_output
    fi
    """
}
