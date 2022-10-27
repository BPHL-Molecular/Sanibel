process readssum {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    #echo ${mypath}
    genome=\$(cat ${pyoutputs} | cut -d "," -f 10)
    #echo \${genome}
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    #echo \${samplename}
    #singularity exec --cleanenv /apps/staphb-toolkit/containers/lyveset_1.1.4f.sif run_assembly_shuffleReads.pl -gz ${mypath}/\${samplename}_1.fq.gz ${mypath}/\${samplename}_2.fq.gz > ${mypath}/\${samplename}_clean_shuffled.fq.gz
    #singularity exec --cleanenv /apps/staphb-toolkit/containers/lyveset_1.1.4f.sif run_assembly_readMetrics.pl ${mypath}/\${samplename}_clean_shuffled.fq.gz -e \${genome} > ${mypath}/\${samplename}_readMetrics.txt
    run_assembly_shuffleReads.pl -gz ${mypath}/\${samplename}_1.fq.gz ${mypath}/\${samplename}_2.fq.gz > ${mypath}/\${samplename}_clean_shuffled.fq.gz
    run_assembly_readMetrics.pl ${mypath}/\${samplename}_clean_shuffled.fq.gz -e \${genome} > ${mypath}/\${samplename}_readMetrics.txt

    """
}
