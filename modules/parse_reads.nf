process parse_reads {
    tag "${meta.id}"

    input:
        tuple val(meta), path(pyoutputs), path(read_metrics)
    output:
        tuple val(meta), path(pyoutputs), emit: out

    script:
    """
    collect_sample_data.py parse-reads \\
        ${read_metrics} \\
        ${pyoutputs}
    """
}
