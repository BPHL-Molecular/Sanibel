process mlst {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/mlst" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("${meta.id}.mlst"), emit: out

    script:
    def prefix = meta.id
    """
    mlst ${assembly} --nopath > ${prefix}.mlst
    """
}
