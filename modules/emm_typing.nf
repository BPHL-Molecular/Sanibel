process emm_typing {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/groupAstrep", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(reads)
    output:
        path("groupAstrep_result.txt"), optional: true
        path("groupAstrep_output/"),    optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})

    if [[ "\${speciesid}" == "Streptococcus pyogenes" || "\${speciesid}" == "Streptococcus dysgalactiae" ]]; then
        emm_typing.py --fastq_1 ${reads[0]} --fastq_2 ${reads[1]} \\
            -m /db/ -o groupAstrep_output > groupAstrep_result.txt
    fi
    """
}
