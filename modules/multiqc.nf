process multiqc {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/multiqc" }, mode: 'copy'

    input:
        tuple val(meta), path(reports)
    output:
        tuple val(meta), path("multiqc_report.html"), emit: report

    script:
    """
    multiqc .
    """
}

process multiqc_global {
    tag "multiqc_global"
    publishDir { "${params.output}" }, mode: 'copy'

    input:
        path summary
        path multiqc_config
        path custom_logo
        path custom_css
        path nf_config
        path mqc_tables, stageAs: 'mqc_in/*'

    output:
        path("sanibel_report.html"), emit: report

    script:
    """
    mkdir -p mqc_in

    {
      echo "# id: 'sanibel_versions'"
      echo "# section_name: 'Software Versions'"
      echo "# description: 'Tool versions from the container tags pinned in nextflow.config.'"
      echo "# plot_type: 'table'"
      echo "# pconfig:"
      echo "#     id: 'sanibel_versions_table'"
      echo "#     col1_header: 'Software'"
      echo "#     no_violin: true"
      echo "#     rows_are_samples: false"
      echo "# headers:"
      echo "#     Version:"
      echo "#         scale: false"
      echo "#         format: '{}'"
      printf 'Software\\tVersion\\n'
      grep -hoE "docker://[^']+" ${nf_config} \\
        | sed -E 's#docker://[^/]*/([^:]+):(.+)#\\1\\t\\2#' \\
        | sort -u
    } > mqc_in/sanibel_versions_mqc.tsv

    multiqc ${params.output} mqc_in \\
        -c ${multiqc_config} \\
        --filename sanibel_report.html \\
        --interactive \\
        --ignore "*/multiqc/*" \\
        --ignore "*/fastqc/*" \\
        --ignore "*sanibel_report*" \\
        --ignore "*sum_report.txt"
    """
}
