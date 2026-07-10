process unicycler {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/assembly" }, mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("${meta.id}.fasta"), emit: assembly

    script:
    def prefix = meta.id
    """
    unicycler \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -o assembly \\
        --min_fasta_length 300 --keep 1 \\
        --min_kmer_frac 0.3 --max_kmer_frac 0.9 \\
        --verbosity 2 --threads ${task.cpus}

    mv assembly/assembly.fasta ${prefix}.fasta
    """
}
