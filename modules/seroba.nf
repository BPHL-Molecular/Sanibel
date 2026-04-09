process seroba {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/seroba", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(reads)
    output:
        path("seroba_output/pred.csv"), optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})
    speciesid2=\$(cut -d "," -f 1 ${pyoutputs})

    if [[ "\${speciesid}" == "Streptococcus pneumoniae" || "\${speciesid2}" == "Streptococcus pneumoniae" ]]; then
        seroba runSerotyping /seroba/database/ ${reads[0]} ${reads[1]} seroba_output
    fi
    """
}
