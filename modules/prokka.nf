process prokka {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/${meta.id}_assembly", mode: 'copy'

    input:
        tuple val(meta), path(assembly), path(pyoutputs)
    output:
        tuple val(meta), path("prokka/"), emit: annotation
        tuple val(meta), path(pyoutputs), path("prokka/${meta.id}.txt"), emit: cds

    script:
    def prefix = meta.id
    """
    genus=\$(cut -d "," -f 1 ${pyoutputs})
    species=\$(cut -d "," -f 2 ${pyoutputs})

    prokka \\
        --genus \${genus} --species \${species} \\
        --strain ${prefix} \\
        --outdir prokka \\
        --prefix ${prefix} \\
        --force --compliant \\
        --locustag \${genus} \\
        --cpus ${task.cpus} \\
        ${assembly}
    """
}
