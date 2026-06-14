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

    # Download one reference genome per candidate in the pool (skip header)
    tail -n +2 ${candidate_pool} | while IFS='\t' read -r species tools accession mash_dist kraken_pct blast16s_pident; do
        safe_name=\$(echo "\${species}" | tr ' ' '_')
        outfile="${prefix}_references/\${safe_name}.fna"

        # Try accession-based download first if accession is known
        if [ "\${accession}" != "NA" ] && [ -n "\${accession}" ]; then
            datasets download genome accession "\${accession}" \\
                --include genome \\
                --filename genome_\${safe_name}.zip 2>/dev/null \\
            && unzip -o genome_\${safe_name}.zip -d genome_\${safe_name}_dir 2>/dev/null \\
            && fna=\$(find genome_\${safe_name}_dir/ncbi_dataset/data -name "*.fna" | head -1) \\
            && [ -n "\${fna}" ] && mv "\${fna}" "\${outfile}" \\
            || true
        fi

        # Fall back to taxon-name search if accession download did not produce a file
        if [ ! -f "\${outfile}" ]; then
            datasets download genome taxon "\${species}" \\
                --reference \\
                --assembly-source refseq \\
                --include genome \\
                --filename genome_\${safe_name}.zip 2>/dev/null \\
            && unzip -o genome_\${safe_name}.zip -d genome_\${safe_name}_dir 2>/dev/null \\
            && fna=\$(find genome_\${safe_name}_dir/ncbi_dataset/data -name "*.fna" | head -1) \\
            && [ -n "\${fna}" ] && mv "\${fna}" "\${outfile}" \\
            || true
        fi

        rm -rf genome_\${safe_name}.zip genome_\${safe_name}_dir 2>/dev/null || true
    done

    # Fail the process only when no references were downloaded at all
    n_downloaded=\$(find ${prefix}_references -name "*.fna" | wc -l)
    [ "\${n_downloaded}" -gt 0 ] || { echo "No reference genomes downloaded for ${prefix}" >&2; exit 1; }
    """
}
