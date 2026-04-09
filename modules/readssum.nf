process readssum {
    tag "${meta.id}"

    input:
        tuple val(meta), path(pyoutputs), path(reads)
    output:
        tuple val(meta), path(pyoutputs), path("${meta.id}_readMetrics.txt"), emit: out

    script:
    def prefix = meta.id
    """
    genome=\$(cut -d "," -f 10 ${pyoutputs})
    run_assembly_shuffleReads.pl -gz ${reads[0]} ${reads[1]} > ${prefix}_clean_shuffled.fq.gz
    run_assembly_readMetrics.pl ${prefix}_clean_shuffled.fq.gz -e \${genome} > ${prefix}_readMetrics.txt
    """
}
