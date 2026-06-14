process download_16s_db {
    storeDir params.blast16s_db_dir

    output:
        path "16s_db", emit: db

    script:
    """
    mkdir -p 16s_db
    cd 16s_db
    wget -q https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
    tar -xzf 16S_ribosomal_RNA.tar.gz
    rm 16S_ribosomal_RNA.tar.gz
    """
}

process blast_16s {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/16S_blast", mode: 'copy'

    input:
        tuple val(meta), path(assembly)
        path db_dir

    output:
        tuple val(meta), path("${meta.id}_16s_blast.tsv"), emit: result

    script:
    def prefix = meta.id
    """
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstitle" \\
        > ${prefix}_16s_blast.tsv
    blastn \\
        -query ${assembly} \\
        -db ${db_dir}/16S_ribosomal_RNA \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
        -max_target_seqs 25 \\
        -max_hsps 1 \\
        -evalue 1e-10 \\
        -perc_identity 80 \\
        -task megablast \\
        -num_threads ${task.cpus} \\
        >> ${prefix}_16s_blast.tsv
    """
}
