process parse_assembly {
    tag "${meta.id}"

    input:
        tuple val(meta), path(distances), path(quast_report)
    output:
        tuple val(meta), path("${meta.id}_pyoutputs.txt"), emit: out

    script:
    def sample = meta.id
    """
    collect_sample_data.py parse-assembly \\
        ${distances} \\
        ${quast_report} \\
        ${sample}_pyoutputs.txt
    """
}
