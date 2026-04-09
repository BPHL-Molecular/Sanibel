process kleborate {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/klebsiella", mode: 'copy'

    input:
        tuple val(meta), path(pyoutputs), path(assembly)
    output:
        path("kleborate_out/klebsiella_pneumo_complex_output.txt"), optional: true
        val meta, emit: done

    script:
    """
    speciesid=\$(cut -d "," -f 21 ${pyoutputs})
    speciesid2=\$(cut -d "," -f 1 ${pyoutputs})

    if [[ "\${speciesid}" == "Klebsiella" || "\${speciesid2}" == "Klebsiella" ]]; then
        kleborate -a ${assembly} -o kleborate_out -p kpsc --trim_headers
    fi
    """
}
