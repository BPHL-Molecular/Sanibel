process mlst {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/mlst" }, mode: 'copy'

    input:
        tuple val(meta), path(assembly)
        path scheme_map
    output:
        tuple val(meta), path("${meta.id}.mlst"), emit: out

    script:
    def prefix  = meta.id
    def species = meta.species ?: 'Unknown'
    def genus   = meta.genus   ?: 'Unknown'
    """
    mlst ${assembly} --nopath > auto.tsv
    scheme=\$(cut -f2 auto.tsv | head -n1)

    tr -d '\\r' < ${scheme_map} > schemes.tsv
    allowed=\$(awk -F'\\t' -v k='${species}' '\$1 == k { print \$2; exit }' schemes.tsv)
    if [[ -z "\$allowed" ]]; then
        allowed=\$(awk -F'\\t' -v k='${genus}' '\$1 == k { print \$2; exit }' schemes.tsv)
    fi

    if [[ '${genus}' == "Unknown" || "\$allowed" == "none" ]]; then
        awk -F'\\t' -v OFS='\\t' '{ print \$1, "-", "-" }' auto.tsv > ${prefix}.mlst
    elif [[ -z "\$allowed" ]]; then
        mv auto.tsv ${prefix}.mlst
    elif [[ ",\$allowed," == *",\$scheme,"* ]]; then
        mv auto.tsv ${prefix}.mlst
    else
        forced=\${allowed%%,*}
        if mlst --list | tr ' ' '\\n' | grep -qx "\$forced"; then
            echo "mlst: auto picked '\$scheme' for ${species}, forcing '\$forced'" >&2
            mlst ${assembly} --scheme "\$forced" --nopath > ${prefix}.mlst
        else
            echo "mlst: WARNING '\$forced' is not in the mlst db, keeping auto '\$scheme'" >&2
            mv auto.tsv ${prefix}.mlst
        fi
    fi
    """
}
