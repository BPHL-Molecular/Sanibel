process kaptive_ab {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/kaptive_ab", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("${meta.id}_ab_k.txt"),  optional: true
        path("${meta.id}_ab_oc.txt"), optional: true
        val meta, emit: done

    script:
    """
    if [[ "${meta.mash_species}" == "Acinetobacter_baumannii" ]]; then
        kaptive.py assembly \\
            /kaptive/reference_database/Acinetobacter_baumannii_k_locus_primary_reference.gbk \\
            ${assembly} -o ${meta.id}_ab_k.txt
        kaptive.py assembly \\
            /kaptive/reference_database/Acinetobacter_baumannii_OC_locus_primary_reference.gbk \\
            ${assembly} -o ${meta.id}_ab_oc.txt
    fi
    """
}
