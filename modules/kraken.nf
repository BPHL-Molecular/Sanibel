process kraken {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/kraken_out", mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("${meta.id}.report"), emit: out

    script:
    def prefix = meta.id
    """
    kraken2 \\
        --db /kraken2-db/minikraken2_v1_8GB/ \\
        --threads ${task.cpus} \\
        --use-names \\
        --report ${prefix}.report \\
        --output ${prefix}_kraken.out \\
        --paired ${reads[0]} ${reads[1]}
    """
}
