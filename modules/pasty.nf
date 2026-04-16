process pasty {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/pasty", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("${meta.id}.tsv"),         optional: true
        path("${meta.id}.details.tsv"), optional: true
        path("${meta.id}.blastn.tsv"),  optional: true
        val meta, emit: done

    script:
    """
    if [[ "${meta.mash_species}" == "Pseudomonas_aeruginosa" ]]; then
        pasty --input ${assembly} --prefix ${meta.id} --outdir .
    fi
    """
}
