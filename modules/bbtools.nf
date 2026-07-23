process bbtools_adapters {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/bbtools" }, mode: 'copy'

    input:
        tuple val(meta), path(trimmed_reads)
    output:
        tuple val(meta), path("${meta.id}_R{1,2}_adaptertrim.fastq.gz"), emit: reads
        tuple val(meta), path("${meta.id}_adapter_stats.txt"),           emit: adapter_stats

    script:
    def prefix = meta.id
    """
    bbduk.sh \\
        in1=${trimmed_reads[0]} \\
        in2=${trimmed_reads[1]} \\
        out1=${prefix}_R1_adaptertrim.fastq.gz \\
        out2=${prefix}_R2_adaptertrim.fastq.gz \\
        ref=/bbmap/resources/adapters.fa \\
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \\
        stats=${prefix}_adapter_stats.txt
    """
}

process bbtools_phix {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/bbtools" }, mode: 'copy'

    input:
        tuple val(meta), path(adaptertrim_reads)
    output:
        tuple val(meta), path("${meta.id}_R{1,2}_clean.fastq.gz"), emit: reads
        tuple val(meta), path("${meta.id}_phix_stats.txt"),         emit: phix_stats

    script:
    def prefix = meta.id
    """
    bbduk.sh \\
        in1=${adaptertrim_reads[0]} \\
        in2=${adaptertrim_reads[1]} \\
        out1=${prefix}_R1_clean.fastq.gz \\
        out2=${prefix}_R2_clean.fastq.gz \\
        ref=/bbmap/resources/phix174_ill.ref.fa.gz \\
        k=31 hdist=1 \\
        stats=${prefix}_phix_stats.txt
    """
}
