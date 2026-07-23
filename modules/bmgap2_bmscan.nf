process bmgap2_bmscan {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy', pattern: 'bmgap2_bmscan'

    input:
        tuple val(meta), path(mlst_file), path(assembly), path(prokka_dir)
    output:
        tuple val(meta), path(mlst_file), emit: out
        tuple val(meta), path("bmgap2_bmscan/${meta.id}_species_analysis.json"), emit: bmscan, optional: true
        path("bmgap2_bmscan"), optional: true

    script:
    """
    mkdir -p assembly_in
    cp -rL ${assembly} ${prokka_dir} assembly_in/

    run_bmgap2_bmscan.py \\
        --sample       ${meta.id} \\
        --mlst         ${mlst_file} \\
        --assembly-dir assembly_in \\
        --outdir       bmgap2_bmscan \\
        --db           ${params.bmgap2_db} \\
        --cpus         ${task.cpus}
    """
}
