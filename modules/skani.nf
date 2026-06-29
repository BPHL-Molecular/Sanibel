process skani {
    tag "${meta.id}"
    publishDir "${params.output}/${meta.id}/skani", mode: 'copy'

    input:
        tuple val(meta), path(assembly), path(references_dir)

    output:
        tuple val(meta), path("${meta.id}_skani.tsv"),         emit: result
        tuple val(meta), path("${meta.id}_skani_species.txt"), emit: species

    script:
    def prefix  = meta.id
    def min_ani = params.skani_routing_min_ani
    """
    skani dist \\
        -t ${task.cpus} \\
        -q ${assembly} \\
        -r ${references_dir}/*.fna \\
        -o ${prefix}_skani.tsv

    # Best-ANI reference species for skani-based routing (Genus_species), emitted only
    # when ANI confirms at the species boundary; empty otherwise so the workflow falls
    # back to Mash. Mirrors parse_skani: top ANI row, ref basename, strip .fna, split __.
    tail -n +2 ${prefix}_skani.tsv | sort -t\$'\\t' -k3 -gr | head -1 \\
        | awk -F'\\t' -v thr=${min_ani} '\$3 >= thr { n=split(\$1,a,"/"); f=a[n]; sub(/\\.fna\$/,"",f); split(f,b,"__"); print b[1] }' \\
        > ${prefix}_skani_species.txt || true
    """
}
