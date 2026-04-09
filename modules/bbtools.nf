process bbtools {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/bbtools", mode: 'copy'

    input:
        tuple val(meta), path(trimmed_reads)
    output:
        tuple val(meta), path("${meta.id}_{1,2}.fq.gz"), emit: reads

    script:
    def prefix = meta.id
    """
    bbduk.sh \\
        in1=${trimmed_reads[0]} in2=${trimmed_reads[1]} \\
        out1=${prefix}_1.rmadpt.fq.gz out2=${prefix}_2.rmadpt.fq.gz \\
        ref=/bbmap/resources/adapters.fa \\
        stats=${prefix}.adapters.stats.txt \\
        ktrim=r k=23 mink=11 hdist=1 tpe tbo \\
        threads=${task.cpus}

    bbduk.sh \\
        in1=${prefix}_1.rmadpt.fq.gz in2=${prefix}_2.rmadpt.fq.gz \\
        out1=${prefix}_1.fq.gz out2=${prefix}_2.fq.gz \\
        outm=${prefix}_matchedphix.fq \\
        ref=/bbmap/resources/phix174_ill.ref.fa.gz \\
        k=31 hdist=1 \\
        stats=${prefix}_phixstats.txt \\
        threads=${task.cpus}
    """
}
