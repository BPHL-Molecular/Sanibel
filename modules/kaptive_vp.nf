process kaptive_vp {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/kaptive_vp", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(assembly)
    output:
        path("${meta.id}_vp_k.txt"), optional: true
        path("${meta.id}_vp_o.txt"), optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})
    speciesid2=\$(cut -d "," -f 1 ${pyoutputs})

    if [[ "\${speciesid}" == "Vibrio parahaemolyticus" || "\${speciesid2}" == "Vibrio parahaemolyticus" ]]; then
        kaptive.py assembly \\
            /kaptive/reference_database/VibrioPara_Kaptivedb_K.gbk \\
            ${assembly} -o ${meta.id}_vp_k.txt
        kaptive.py assembly \\
            /kaptive/reference_database/VibrioPara_Kaptivedb_O.gbk \\
            ${assembly} -o ${meta.id}_vp_o.txt
    fi
    """
}
