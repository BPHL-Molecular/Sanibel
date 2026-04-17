process kaptive_vp {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/kaptive_vp", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("${meta.id}_vp_k.txt")
        path("${meta.id}_vp_o.txt")
        val meta, emit: done

    script:
    """
    kaptive.py assembly \\
        /kaptive/reference_database/VibrioPara_Kaptivedb_K.gbk \\
        ${assembly} -o ${meta.id}_vp_k.txt
    kaptive.py assembly \\
        /kaptive/reference_database/VibrioPara_Kaptivedb_O.gbk \\
        ${assembly} -o ${meta.id}_vp_o.txt
    """
}
