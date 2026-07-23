process seqsero2 {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/seqsero2" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        path("SeqSero_result_*")
        val meta, emit: done

    script:
    """
    SeqSero2_package.py -p ${task.cpus} -t 2 -i ${reads[0]} ${reads[1]}
    """
}
