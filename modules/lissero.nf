process lissero {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/lissero" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("lissero_output.txt")
        val meta, emit: done

    script:
    """
    lissero ${assembly} > lissero_output.txt
    """
}
