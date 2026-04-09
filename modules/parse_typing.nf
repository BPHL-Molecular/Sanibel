// thin wrapper — logic lives in bin/collect_sample_data.py parse-typing
process parse_typing {
    tag "${meta.id}"

    input:
        tuple val(meta), path(pyoutputs), path(cds_txt), path(mlst_out), path(kraken_report), path(pmga_txt), path(neisseria_txt), path(hinfluenzae_txt)
    output:
        tuple val(meta), path(pyoutputs), emit: out

    script:
    def sample = meta.id
    def mypath = "${params.output}/${meta.id}"
    """
    collect_sample_data.py parse-typing \\
        --prokka_txt ${cds_txt} \\
        --mlst_file ${mlst_out} \\
        --kraken_report ${kraken_report} \\
        --pmga_txt ${pmga_txt} \\
        --neisseria_txt ${neisseria_txt} \\
        --hinfluenzae_txt ${hinfluenzae_txt} \\
        --mypath ${mypath} \\
        --sample_id ${sample} \\
        ${pyoutputs}
    """
}
