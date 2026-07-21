process bmgap2_amr {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy', pattern: 'bmgap2_amr'

    input:
        tuple val(meta), path(mlst_file), path(pmga_files)
    output:
        tuple val(meta), path(mlst_file), emit: out
        tuple val(meta), path("bmgap2_amr/${meta.id}*amr_data.json"), emit: amr, optional: true
        path("bmgap2_amr"), optional: true

    script:
    """
    run_bmgap2_amr.py \\
        --sample   ${meta.id} \\
        --mlst     ${mlst_file} \\
        --pmga-dir . \\
        --outdir   bmgap2_amr \\
        --db       ${params.bmgap2_db}
    """
}
