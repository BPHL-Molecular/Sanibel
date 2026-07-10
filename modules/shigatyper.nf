process shigatyper {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/shigatyper" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        path("shigatyper_output.txt")
        val meta, emit: done

    script:
    """
    shigatyper --R1 ${reads[0]} --R2 ${reads[1]} > shigatyper_output.txt
    """
}
