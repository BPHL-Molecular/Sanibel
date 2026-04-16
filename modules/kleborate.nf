process kleborate {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/kleborate", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
    output:
        path("kleborate_out/klebsiella_pneumo_complex_output.txt"), optional: true
        val meta, emit: done

    script:
    """
    if [[ "${meta.mash_genus}" == "Klebsiella" ]]; then
        kleborate -a ${assembly} -o kleborate_out -p kpsc --trim_headers
    fi
    """
}
