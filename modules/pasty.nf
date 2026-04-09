process pasty {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/pasty", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(assembly)
    output:
        path("${meta.id}.tsv"),         optional: true
        path("${meta.id}.details.tsv"), optional: true
        path("${meta.id}.blastn.tsv"),  optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})
    speciesid2=\$(cut -d "," -f 1 ${pyoutputs})

    if [[ "\${speciesid}" == "Pseudomonas aeruginosa" || "\${speciesid2}" == "Pseudomonas aeruginosa" ]]; then
        pasty --input ${assembly} --prefix ${meta.id} --outdir .
    fi
    """
}
