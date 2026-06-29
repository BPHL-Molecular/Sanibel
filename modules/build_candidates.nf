process build_candidates {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/candidate_pool", mode: 'copy'

    input:
        tuple val(meta), path(distances), path(kraken_report), path(blast16s)

    output:
        tuple val(meta), path("${meta.id}_candidate_pool.tsv"), emit: pool

    script:
    def prefix = meta.id
    """
    build_candidate_pool.py \\
        ${distances} \\
        ${kraken_report} \\
        ${blast16s} \\
        > ${prefix}_candidate_pool.tsv
    """
}
