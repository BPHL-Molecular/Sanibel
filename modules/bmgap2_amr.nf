process bmgap2_amr {
    tag "${meta.id}"

    input:
        tuple val(meta), path(pyoutputs)
    output:
        tuple val(meta), path(pyoutputs), emit: out

    script:
    def mypath = "${params.output}/${meta.id}"
    """
    run_bmgap2_amr.py ${mypath} ${pyoutputs} ${params.bmgap2_db}
    """
}

