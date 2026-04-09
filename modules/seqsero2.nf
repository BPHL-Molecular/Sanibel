process seqsero2 {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/salmonella", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(reads)
    output:
        path("SeqSero_result_*"), optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})
    speciesid2=\$(cut -d "," -f 1 ${pyoutputs})

    if [[ "\${speciesid}" == "Salmonella" || "\${speciesid2}" == "Salmonella" ]]; then
        SeqSero2_package.py -p ${task.cpus} -t 2 -i ${reads[0]} ${reads[1]}
    fi
    """
}
