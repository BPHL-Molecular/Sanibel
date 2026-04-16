process parse_assembly {
    tag "${meta.id}"

    input:
        tuple val(meta), path(distances), path(quast_report)
    output:
        tuple val(meta), path("${meta.id}_assembly_stats.txt"), emit: out

    script:
    def sample = meta.id
    """
    collect_sample_data.py parse-assembly \\
        ${distances} \\
        ${quast_report} \\
        ${sample}_assembly_stats.txt
    """
}
