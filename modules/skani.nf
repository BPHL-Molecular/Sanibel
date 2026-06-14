process skani {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/skani_ANI", mode: 'copy'

    input:
        tuple val(meta), path(assembly), path(references_dir)

    output:
        tuple val(meta), path("${meta.id}_skani.tsv"), emit: result

    script:
    def prefix = meta.id
    """
    skani dist \\
        -t ${task.cpus} \\
        -q ${assembly} \\
        -r ${references_dir}/*.fna \\
        -o ${prefix}_skani.tsv
    """
}
