process legsta {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/legsta" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("legsta_output.txt")
        val meta, emit: done

    script:
    """
    legsta ${assembly} > legsta_output.txt
    """
}
