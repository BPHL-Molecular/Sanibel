process fastqc2 {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${x}"
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   
   
   fastqc ${params.output}/${x}/${x}_1.fq.gz ${params.output}/${x}/${x}_2.fq.gz

   mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_clean_fastqc.html
   mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_clean_fastqc.zip
   mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_clean_fastqc.html
   mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_clean_fastqc.zip
   
   

   """
}
