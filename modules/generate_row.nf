process generate_row {
    tag "${meta.id}"

    input:
        tuple val(meta), path(pyoutputs), path(hinfluenzae_txt)
    output:
        tuple val(meta), path("${meta.id}_row*.tsv"), emit: row

    script:
    def sample = meta.id
    def mypath = "${params.output}/${meta.id}"
    """
    generate_report_row.py ${pyoutputs} ${sample} ${mypath} \\
        --meningitis ${params.meningitis} \\
        --hinfluenzae_txt ${hinfluenzae_txt}
    """
}

