process kaptive {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/kaptive_${variant}", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
        val variant
    output:
        path("${meta.id}_${variant}_*.txt")
        val meta, emit: done

    script:
    // Per-variant locus/serotype databases and output names. Output filenames are
    // unchanged from the original kaptive_ab / kaptive_vp modules so summary_report.py
    // still finds them (note ab uses the _oc suffix, vp uses _o).
    def cfg = [
        ab: [ k: [ db: 'ab_k', out: "${meta.id}_ab_k.txt"  ],
              o: [ db: 'ab_o', out: "${meta.id}_ab_oc.txt" ] ],
        vp: [ k: [ db: '/kaptive/reference_database/VibrioPara_Kaptivedb_K.gbk', out: "${meta.id}_vp_k.txt" ],
              o: [ db: '/kaptive/reference_database/VibrioPara_Kaptivedb_O.gbk', out: "${meta.id}_vp_o.txt" ] ],
    ][variant]
    """
    kaptive assembly ${cfg.k.db} ${assembly} -o ${cfg.k.out}
    kaptive assembly ${cfg.o.db} ${assembly} -o ${cfg.o.out}
    """
}
