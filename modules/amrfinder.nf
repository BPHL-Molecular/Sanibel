process amrfinder {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/${meta.id}_assembly", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("${meta.id}_amrfinderplus_report.tsv"), emit: out

    script:
    def prefix = meta.id
    """
    amrfinder \\
        --num_threads ${task.cpus} \\
        -n ${assembly} \\
        --plus \\
        -o ${prefix}_amrfinderplus_report.tsv
    """
}
