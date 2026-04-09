process plasmidfinder {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/plasmid", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(reads)
    output:
        path("results_tab.tsv"), optional: true
        val meta, emit: done

    script:
    """
    python -m plasmidfinder -i ${reads[0]} ${reads[1]} -o ./
    """
}
