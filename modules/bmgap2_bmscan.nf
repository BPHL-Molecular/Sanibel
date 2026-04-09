process bmgap2_bmscan {
    tag "${meta.id}"

    input:
        tuple val(meta), path(pyoutputs)
    output:
        tuple val(meta), path(pyoutputs), emit: out

    script:
    def mypath = "${params.output}/${meta.id}"
    """
    run_bmgap2_bmscan.py ${mypath} ${pyoutputs} ${params.bmgap2_db} ${task.cpus}
    """
}

