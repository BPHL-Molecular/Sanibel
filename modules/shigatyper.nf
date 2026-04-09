process shigatyper {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/shigatyper", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(reads)
    output:
        path("shigatyper_output.txt"), optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})
    speciesid2=\$(cut -d "," -f 1 ${pyoutputs})

    if [[ "\${speciesid}" == "Shigella" || "\${speciesid2}" == "Shigella" ]]; then
        shigatyper --R1 ${reads[0]} --R2 ${reads[1]} > shigatyper_output.txt
    fi
    """
}
