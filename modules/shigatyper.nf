process shigatyper {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/shigatyper", mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        path("shigatyper_output.txt"), optional: true
        val meta, emit: done

    script:
    """
    if [[ "${meta.mash_genus}" == "Shigella" ]]; then
        shigatyper --R1 ${reads[0]} --R2 ${reads[1]} > shigatyper_output.txt
    fi
    """
}
