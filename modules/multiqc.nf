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
        path nf_config

    output:
        path("interactive_report.html"), emit: report

    script:
    """
    # Tool versions from the pinned container tags in nextflow.config
    {
      echo "software_versions:"
      grep -hoE "docker://[^']+" ${nf_config} \\
        | sed -E 's#docker://[^/]*/([^:]+):(.+)#  \\1: "\\2"#; s/^  kraken2:/  kraken:/; s/^  bbtools:/  bbmap:/' \\
        | sort -u
    } > sanibel_versions.yml

    multiqc ${params.output} \\
        -c ${multiqc_config} \\
        -c sanibel_versions.yml \\
        --filename interactive_report.html \\
        --interactive \\
        --ignore "*/multiqc/*" \\
        --ignore "*interactive_report*" \\
        --ignore "*sum_report.txt"
    """
}
