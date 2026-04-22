process kaptive_ab {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/kaptive_ab", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("${meta.id}_ab_k.txt")
        path("${meta.id}_ab_oc.txt")
        val meta, emit: done

    script:
    """
    kaptive assembly ab_k ${assembly} -o ${meta.id}_ab_k.txt
    kaptive assembly ab_o ${assembly} -o ${meta.id}_ab_oc.txt
    """
}
