process mash {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/mash" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("${meta.id}_mash_distances.tab"), emit: distances

    script:
    def prefix = meta.id
    """
    zcat ${reads[0]} ${reads[1]} | mash sketch -p ${task.cpus} -r -m 2 - -o ${prefix}_sketch

    mash dist -p ${task.cpus} /db/RefSeqProkSketchesv235.msh ${prefix}_sketch.msh > ${prefix}_distances.tab
    sort -gk3 ${prefix}_distances.tab | head -50 > ${prefix}_mash_distances.tab
    """
}
