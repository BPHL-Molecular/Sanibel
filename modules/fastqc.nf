process fastqc {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   mkdir -p ${params.output}/${x}
   cp ${params.input}/${x}_*.fastq.gz ${params.output}/${x}

   ##### pull container from docker hub, and if multiple containers in a process
   fastqc ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz

   mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_original_fastqc.html
   mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_original_fastqc.zip
   mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_original_fastqc.html
   mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_original_fastqc.zip
   
   """
}
