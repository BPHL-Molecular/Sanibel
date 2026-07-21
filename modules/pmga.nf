process pmga {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/pmga" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly), path(mlst_out)
    output:
        tuple val(meta), path("${meta.id}sta.txt"), emit: out
        tuple val(meta), path("${meta.id}sta*"),    emit: files

    script:
    def prefix = meta.id
    """
    scheme=\$(awk 'NR==1{print \$2}' ${mlst_out})
    pmga --blastdir /pmga/blastdbs -o pmga_out --force \\
        --threads ${task.cpus} \\
        --species \${scheme} ${assembly}
    cp pmga_out/${prefix}sta* .
    """
}
