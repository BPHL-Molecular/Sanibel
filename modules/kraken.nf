process kraken {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    
    mkdir -p ${mypath}/kraken_out/
    #singularity exec --cleanenv /apps/staphb-toolkit/containers/kraken2_2.0.8-beta.sif kraken2 --db /kraken2-db/minikraken2_v1_8GB/ --use-names --report ${mypath}/kraken_out/\${samplename}.report --output ${mypath}/kraken_out/\${samplename}_kraken.out --paired ${params.input}/\${samplename}_1.fastq.gz ${params.input}/\${samplename}_2.fastq.gz
    kraken2 --db /kraken2-db/minikraken2_v1_8GB/ --use-names --report ${mypath}/kraken_out/\${samplename}.report --output ${mypath}/kraken_out/\${samplename}_kraken.out --paired ${params.input}/\${samplename}_1.fastq.gz ${params.input}/\${samplename}_2.fastq.gz

    mkdir -p assemblies
    mkdir -p annotations
    mkdir -p amrfinder_results
    cp ${mypath}/\${samplename}_assembly/\${samplename}.fasta assemblies/
    cp ${mypath}/\${samplename}_assembly/prokka/\${samplename}.gff annotations/
    cp ${mypath}/\${samplename}_assembly/\${samplename}_amrfinderplus_report.tsv amrfinder_results/ 
    """
}
