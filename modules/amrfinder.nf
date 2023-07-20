process amrfinder {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    
    #singularity exec --cleanenv -B /tmp:/tmp -B ${mypath}/\${samplename}_assembly/:/data docker://staphb/ncbi-amrfinderplus:3.10.1 amrfinder -n /data/\${samplename}.fasta --plus -o /data/\${samplename}_amrfinderplus_report.tsv
    amrfinder -n ${mypath}/\${samplename}_assembly/\${samplename}.fasta --plus -o ${mypath}/\${samplename}_assembly/\${samplename}_amrfinderplus_report.tsv

    """
}
