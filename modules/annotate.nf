process annotate {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    genus=\$(cat ${pyoutputs} | cut -d "," -f 1)
    species=\$(cat ${pyoutputs} | cut -d "," -f 2)
    echo \${genus}
    #singularity exec --cleanenv /apps/staphb-toolkit/containers/prokka_1.14.5.sif prokka --genus \${genus} --species \${species} --strain \${samplename} --outdir ${mypath}/\${samplename}_assembly/prokka --prefix \${samplename} --force --compliant --locustag \${genus} ${mypath}/\${samplename}_assembly/\${samplename}.fasta
    singularity exec --cleanenv docker://staphb/prokka:1.14.5 prokka --genus \${genus} --species \${species} --strain \${samplename} --outdir ${mypath}/\${samplename}_assembly/prokka --prefix \${samplename} --force --compliant --locustag \${genus} ${mypath}/\${samplename}_assembly/\${samplename}.fasta

    #singularity exec --cleanenv -B /tmp:/tmp -B ${mypath}/\${samplename}_assembly/:/data /apps/staphb-toolkit/containers/ncbi-amrfinderplus_3.10.1.sif amrfinder -n /data/\${samplename}.fasta --plus -o /data/\${samplename}_amrfinderplus_report.tsv
    singularity exec --cleanenv -B /tmp:/tmp -B ${mypath}/\${samplename}_assembly/:/data docker://staphb/ncbi-amrfinderplus:3.10.1 amrfinder -n /data/\${samplename}.fasta --plus -o /data/\${samplename}_amrfinderplus_report.tsv
    
    #singularity exec --bind ${mypath}/\${samplename}_assembly:/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/mlst_2.19.0.sif mlst ${mypath}/\${samplename}_assembly/\${samplename}.fasta --nopath > ${mypath}/\${samplename}_assembly/\${samplename}.mlst
    singularity exec --bind ${mypath}/\${samplename}_assembly:/data --pwd /data --cleanenv docker://staphb/mlst:2.19.0 mlst ${mypath}/\${samplename}_assembly/\${samplename}.fasta --nopath > ${mypath}/\${samplename}_assembly/\${samplename}.mlst

    mkdir -p ${mypath}/kraken_out/
    #singularity exec --cleanenv /apps/staphb-toolkit/containers/kraken2_2.0.8-beta.sif kraken2 --db /kraken2-db/minikraken2_v1_8GB/ --use-names --report ${mypath}/kraken_out/\${samplename}.report --output ${mypath}/kraken_out/\${samplename}_kraken.out --paired ${params.input}/\${samplename}_1.fastq.gz ${params.input}/\${samplename}_2.fastq.gz
    singularity exec --cleanenv docker://staphb/kraken2:2.0.8-beta kraken2 --db /kraken2-db/minikraken2_v1_8GB/ --use-names --report ${mypath}/kraken_out/\${samplename}.report --output ${mypath}/kraken_out/\${samplename}_kraken.out --paired ${params.input}/\${samplename}_1.fastq.gz ${params.input}/\${samplename}_2.fastq.gz

    mkdir -p assemblies
    mkdir -p annotations
    mkdir -p amrfinder_results
    cp ${mypath}/\${samplename}_assembly/\${samplename}.fasta assemblies/
    cp ${mypath}/\${samplename}_assembly/prokka/\${samplename}.gff annotations/
    cp ${mypath}/\${samplename}_assembly/\${samplename}_amrfinderplus_report.tsv amrfinder_results/ 
    """
}
