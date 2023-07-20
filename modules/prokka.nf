process prokka {
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
    
    prokka --genus \${genus} --species \${species} --strain \${samplename} --outdir ${mypath}/\${samplename}_assembly/prokka --prefix \${samplename} --force --compliant --locustag \${genus} ${mypath}/\${samplename}_assembly/\${samplename}.fasta

     
    """
}
