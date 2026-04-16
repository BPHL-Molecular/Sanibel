process legsta {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/legsta", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("legsta_output.txt"), optional: true
        val meta, emit: done

    script:
    """
    if [[ "${meta.mash_species}" == "Legionella_pneumophila" ]]; then
        legsta ${assembly} > legsta_output.txt
    fi
    """
}
