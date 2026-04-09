process serotypefinder {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/escherichia", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(reads)
    output:
        path("results_tab.tsv"), optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})

    if [[ "\${speciesid}" == "Escherichia coli" ]]; then
        serotypefinder.py -i ${reads[0]} ${reads[1]} -o ./
    fi
    """
}
