process prokka {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/assembly", mode: 'copy'

    input:
        tuple val(meta), path(assembly), path(assembly_stats)
    output:
        tuple val(meta), path("prokka/"), emit: annotation
        tuple val(meta), path("prokka/${meta.id}.txt"), emit: cds_txt

    script:
    def prefix = meta.id
    """
    export _JAVA_OPTIONS="-XX:-UsePerfData"
    genus=\$(cut -d "," -f 1 ${assembly_stats})
    species=\$(cut -d "," -f 2 ${assembly_stats})

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
