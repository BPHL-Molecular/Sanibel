process prokka {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/assembly" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("prokka/"), emit: annotation
        tuple val(meta), path("prokka/${meta.id}.txt"), emit: cds_txt

    script:
    def prefix  = meta.id
    // skani-confirmed call, falling back to Mash when skani could not resolve one
    def taxon   = (meta.species && meta.species != 'Unknown') ? meta.species
                                                              : (meta.mash_species ?: 'Unknown_unknown')
    def parts   = taxon.tokenize('_')
    def genus   = parts[0]
    def species = parts.size() > 1 ? parts[1] : 'unknown'
    """
    export _JAVA_OPTIONS="-XX:-UsePerfData"

    prokka \\
        --genus ${genus} --species ${species} \\
        --strain ${prefix} \\
        --outdir prokka \\
        --prefix ${prefix} \\
        --force --compliant \\
        --locustag ${genus} \\
        --cpus ${task.cpus} \\
        ${assembly}
    """
}
