process mlst {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    
    #singularity exec --bind ${mypath}/\${samplename}_assembly:/data --pwd /data --cleanenv docker://staphb/mlst:2.19.0 mlst ${mypath}/\${samplename}_assembly/\${samplename}.fasta --nopath > ${mypath}/\${samplename}_assembly/\${samplename}.mlst
    mlst ${mypath}/\${samplename}_assembly/\${samplename}.fasta --nopath > ${mypath}/\${samplename}_assembly/\${samplename}.mlst

    
    """
}
