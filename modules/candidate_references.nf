process candidate_references {
    tag "${meta.id}"
    publishDir { "${params.output}/${meta.id}/reference_genomes" }, mode: 'copy', pattern: '*_reference_genomes.txt'

    input:
        tuple val(meta), path(candidate_pool)

    output:
        tuple val(meta), path("${meta.id}_references/"),      emit: references
        path("${meta.id}_reference_genomes.txt"),             emit: manifest

    script:
    def prefix = meta.id
    def n_refs = 5
    """
    mkdir -p ${prefix}_references
    export REF_DIR="\$PWD/${prefix}_references"
    export N_REFS="${n_refs}"

    extract_all() {
        unz_dir="\$1"; safe_name="\$2"
        for accdir in "\${unz_dir}/ncbi_dataset/data"/*/; do
            [ -d "\${accdir}" ] || continue
            acc=\$(basename "\${accdir}")
            fna=\$(find "\${accdir}" -name "*.fna" | head -1)
            [ -n "\${fna}" ] && cp "\${fna}" "\${REF_DIR}/\${safe_name}__\${acc}.fna"
        done
    }
    export -f extract_all

    retry_ncbi() {
        attempt=1
        while true; do
            "\$@" && return 0
            [ \${attempt} -ge 4 ] && return 1
            sleep \$(( attempt * attempt * 5 ))
            attempt=\$(( attempt + 1 ))
        done
    }
    export -f retry_ncbi

    download_candidate() {
        species="\$1"
        accession="\$2"
        safe_name=\$(echo "\${species}" | tr ' ' '_')
        workdir="dl_\${safe_name}"

        # 1. Non-zero exit is an NCBI failure; clean exit with no accessions means the species has none
        if retry_ncbi datasets summary genome taxon "\${species}" \\
                --assembly-source refseq \\
                --limit "\${N_REFS}" \\
                --report ids_only --as-json-lines \\
                > "sum_\${safe_name}.jsonl" 2> "sum_\${safe_name}.err"; then
            acc_list=\$(grep -oE 'GC[FA]_[0-9]+\\.[0-9]+' "sum_\${safe_name}.jsonl" || true)
        else
            echo "candidate_references: datasets summary failed for \${species}" >&2
            cat "sum_\${safe_name}.err" >&2
            echo "\${species}" >> ncbi_failures.txt
            return 0
        fi

        # 2. Always include the sketch representative
        if [ "\${accession}" != "NA" ] && [ -n "\${accession}" ]; then
            acc_list=\$(printf '%s\\n%s\\n' "\${acc_list}" "\${accession}")
        fi
        acc_list=\$(printf '%s\\n' "\${acc_list}" | grep -oE 'GC[FA]_[0-9]+\\.[0-9]+' | sort -u || true)

        # 3. Download the selected accessions in a single call and split per accession
        if [ -n "\${acc_list}" ]; then
            if retry_ncbi datasets download genome accession \${acc_list} \\
                    --include genome \\
                    --filename "\${workdir}.zip" 2> "dl_\${safe_name}.err"; then
                unzip -o "\${workdir}.zip" -d "\${workdir}_dir" >/dev/null 2>&1 \\
                    && extract_all "\${workdir}_dir" "\${safe_name}" || true
            else
                echo "candidate_references: datasets download failed for \${species}" >&2
                cat "dl_\${safe_name}.err" >&2
                echo "\${species}" >> ncbi_failures.txt
            fi
            rm -rf "\${workdir}.zip" "\${workdir}_dir" 2>/dev/null || true
        fi

        # 4. Last-resort fallback (summary/download produced nothing): old reference path
        if ! ls "\${REF_DIR}/\${safe_name}__"*.fna >/dev/null 2>&1; then
            echo "candidate_references: reference fallback for \${species}, expect one genome only" >&2
            refwork="ref_\${safe_name}"
            if retry_ncbi datasets download genome taxon "\${species}" \\
                    --reference \\
                    --assembly-source refseq \\
                    --include genome \\
                    --filename "\${refwork}.zip" 2> "ref_\${safe_name}.err"; then
                unzip -o "\${refwork}.zip" -d "\${refwork}_dir" >/dev/null 2>&1 \\
                    && extract_all "\${refwork}_dir" "\${safe_name}" || true
            else
                echo "candidate_references: reference fallback failed for \${species}" >&2
                cat "ref_\${safe_name}.err" >&2
                echo "\${species}" >> ncbi_failures.txt
            fi
            rm -rf "\${refwork}.zip" "\${refwork}_dir" 2>/dev/null || true
        fi
    }
    export -f download_candidate

    tail -n +2 ${candidate_pool} | while IFS=\$'\\t' read -r species tools accession rest; do
        [ -n "\${species}" ] && printf '%s\\0%s\\0' "\${species}" "\${accession}"
    done | xargs -0 -n 2 -P 2 bash -c 'download_candidate "\$0" "\$1"'

    # Manifest of the references skani will use (record only; the genomes themselves are not published)
    ls ${prefix}_references/*.fna 2>/dev/null | xargs -r -n1 basename \\
        | sed 's/\\.fna\$//; s/__/_/' | sort > ${prefix}_reference_genomes.txt

    if [ -s ncbi_failures.txt ]; then
        sort -u ncbi_failures.txt > failed_candidates.txt
        echo "candidate_references: NCBI lookups failed for: \$(tr '\\n' ' ' < failed_candidates.txt)" >&2
        tail -n +2 ${candidate_pool} | head -3 | cut -f1 > top_candidates.txt
        if grep -qxF -f top_candidates.txt failed_candidates.txt; then
            echo "candidate_references: a top-ranked candidate failed, failing the task so it retries" >&2
            exit 1
        fi
    fi

    # Fail the process only when no references were downloaded at all
    n_downloaded=\$(find ${prefix}_references -name "*.fna" | wc -l)
    [ "\${n_downloaded}" -gt 0 ] || { echo "No reference genomes downloaded for ${prefix}" >&2; exit 1; }
    """
}
