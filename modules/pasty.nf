process pasty {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/pasty", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("${meta.id}.tsv")
        path("${meta.id}.details.tsv")
        path("${meta.id}.blastn.tsv")
        val meta, emit: done

    script:
    """
    pasty --input ${assembly} --prefix ${meta.id} --outdir .
    """
}
