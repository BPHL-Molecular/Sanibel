process bmgap2_bmscan {
    tag "${meta.id}"

    input:
        tuple val(meta), path(mlst_file)
    output:
        tuple val(meta), path(mlst_file), emit: out

    script:
    def mypath = "${params.output}/${meta.id}"
    """
    run_bmgap2_bmscan.py ${mypath} ${mlst_file} ${params.bmgap2_db} ${task.cpus}
    """
}

