process bmgap2_locusextractor {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}" }, mode: 'copy', pattern: 'bmgap2_locusextractor'

    input:
        tuple val(meta), path(mlst_file), path(assembly), path(prokka_dir)
    output:
        tuple val(meta), path(mlst_file), emit: out
        tuple val(meta), path("bmgap2_locusextractor/LE_*"), emit: le, optional: true
        path("bmgap2_locusextractor"), optional: true

    script:
    """
    mkdir -p assembly_in
    cp -rL ${assembly} ${prokka_dir} assembly_in/

    run_bmgap2_locusextractor.py \\
        --sample       ${meta.id} \\
        --mlst         ${mlst_file} \\
        --assembly-dir assembly_in \\
        --outdir       bmgap2_locusextractor \\
        --db           ${params.bmgap2_db}
    """
}
