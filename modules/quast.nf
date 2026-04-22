process quast {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/assembly", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        tuple val(meta), path("quast_results/report.tsv"), emit: report

    script:
    """
    quast.py --threads ${task.cpus} -o quast_results ${assembly}
    """
}
