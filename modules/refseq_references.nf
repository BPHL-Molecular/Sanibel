process refseq_references {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/refseq_references", mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple val(meta), path(candidate_pool)

    output:
        tuple val(meta), path("${meta.id}_references/"),      emit: references

    script:
    def prefix = meta.id
    """
    mkdir -p ${prefix}_references
    export REF_DIR="\$PWD/${prefix}_references"

    # Download one reference genome for a single candidate: accession first,
    # then taxon-name fallback. Failures are tolerated (handled by the final guard).
    download_one() {
        species="\$1"
        accession="\$2"
        safe_name=\$(echo "\${species}" | tr ' ' '_')
        outfile="\${REF_DIR}/\${safe_name}.fna"
        workdir="dl_\${safe_name}"

        # Try accession-based download first if accession is known
        if [ "\${accession}" != "NA" ] && [ -n "\${accession}" ]; then
            datasets download genome accession "\${accession}" \\
                --include genome \\
                --filename "\${workdir}.zip" 2>/dev/null \\
            && unzip -o "\${workdir}.zip" -d "\${workdir}_dir" 2>/dev/null \\
            && fna=\$(find "\${workdir}_dir/ncbi_dataset/data" -name "*.fna" | head -1) \\
            && [ -n "\${fna}" ] && mv "\${fna}" "\${outfile}" \\
            || true
        fi

        # Fall back to taxon-name search if accession download did not produce a file
        if [ ! -f "\${outfile}" ]; then
            datasets download genome taxon "\${species}" \\
                --reference \\
                --assembly-source refseq \\
                --include genome \\
                --filename "\${workdir}.zip" 2>/dev/null \\
            && unzip -o "\${workdir}.zip" -d "\${workdir}_dir" 2>/dev/null \\
            && fna=\$(find "\${workdir}_dir/ncbi_dataset/data" -name "*.fna" | head -1) \\
            && [ -n "\${fna}" ] && mv "\${fna}" "\${outfile}" \\
            || true
        fi

        rm -rf "\${workdir}.zip" "\${workdir}_dir" 2>/dev/null || true
    }
    export -f download_one

    # Feed each candidate (skip header) as a NUL-delimited species/accession pair
    # to xargs, running up to 3 downloads in parallel.
    tail -n +2 ${candidate_pool} | while IFS=\$'\\t' read -r species tools accession rest; do
        [ -n "\${species}" ] && printf '%s\\0%s\\0' "\${species}" "\${accession}"
    done | xargs -0 -n 2 -P 3 bash -c 'download_one "\$0" "\$1"'

    # Fail the process only when no references were downloaded at all
    n_downloaded=\$(find ${prefix}_references -name "*.fna" | wc -l)
    [ "\${n_downloaded}" -gt 0 ] || { echo "No reference genomes downloaded for ${prefix}" >&2; exit 1; }
    """
}
