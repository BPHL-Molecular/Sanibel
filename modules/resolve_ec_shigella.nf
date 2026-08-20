process resolve_ec_shigella {
    tag "${meta.id}"

    input:
        tuple val(meta), path(skani_tsv), path(shigatyper_txt)

    output:
        tuple val(meta), path("${meta.id}_species_resolved.txt"), emit: resolved

    script:
    """
    resolve_ec_shigella.py \\
        --skani      ${skani_tsv} \\
        --shigatyper ${shigatyper_txt} \\
        --prefix     ${meta.id}
    """
}
