process multiqc {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${x}"
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   
   multiqc ${params.output}/${x}/${x}_*_fastqc.zip -o ${params.output}/${x}
     
   """
}
