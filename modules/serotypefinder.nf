process serotypefinder {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/serotypefinder", mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        path("results_tab.tsv")
        val meta, emit: done

    script:
    """
    serotypefinder.py -i ${reads[0]} ${reads[1]} -o ./
    """
}
