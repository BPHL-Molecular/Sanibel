process shigatyper {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/shigatyper" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("shigatyper_output.txt"), emit: result
        val meta, emit: done

    script:
    """
    shigatyper --R1 ${reads[0]} --R2 ${reads[1]} > shigatyper_output.txt
    """
}
