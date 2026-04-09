process legsta {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/legsta", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(assembly)
    output:
        path("legsta_output.txt"), optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})
    speciesid2=\$(cut -d "," -f 1 ${pyoutputs})

    if [[ "\${speciesid}" == "Legionella pneumophila" || "\${speciesid2}" == "Legionella pneumophila" ]]; then
        legsta ${assembly} > legsta_output.txt
    fi
    """
}
