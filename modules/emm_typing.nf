process emm_typing {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/emm_typing", mode: 'copy'
    containerOptions "--bind ${task.workDir}:/EMBOSS-6.6.0/emboss/.libs"

    input:
        tuple val(meta), path(reads)
    output:
        path("groupAstrep_result.txt")
        path("groupAstrep_output/")
        val meta, emit: done

    script:
    """
    mkdir -p groupAstrep_output
    emm_typing.py --fastq_1 ${reads[0]} --fastq_2 ${reads[1]} \\
        -m /db/ -o groupAstrep_output > groupAstrep_result.txt
    """
}
