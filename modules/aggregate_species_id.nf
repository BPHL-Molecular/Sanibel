process aggregate_species_id {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/candidate_species", mode: 'copy'

    input:
        tuple val(meta), path(assembly_stats), path(kraken_report), path(blast_16s_tsv)

    output:
        tuple val(meta), path("${meta.id}_candidate_species.txt"), emit: out

    script:
    def prefix = meta.id
    """
    aggregate_species_id.py \\
        ${assembly_stats} \\
        ${kraken_report} \\
        ${blast_16s_tsv} \\
        > ${prefix}_candidate_species.txt
    """
}
