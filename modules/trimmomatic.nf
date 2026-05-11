process trimmomatic {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/trimmomatic", mode: 'copy'

    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("${meta.id}_trim_{1,2}.fastq.gz"), emit: reads

    script:
    def prefix = meta.id
    """
    trimmomatic PE -threads ${task.cpus} -phred33 \\
        -trimlog ${prefix}.log \\
        ${reads[0]} ${reads[1]} \\
        ${prefix}_trim_1.fastq.gz ${prefix}_unpaired_trim_1.fastq.gz \\
        ${prefix}_trim_2.fastq.gz ${prefix}_unpaired_trim_2.fastq.gz \\
        ILLUMINACLIP:/usr/local/bin/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:71 TRAILING:20
    """
}
