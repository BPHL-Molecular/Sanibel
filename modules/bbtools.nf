process bbtools {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${x}"
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   
   bbduk.sh in1=${params.output}/${x}/${x}_trim_1.fastq.gz in2=${params.output}/${x}/${x}_trim_2.fastq.gz out1=${params.output}/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo

   bbduk.sh in1=${params.output}/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/${x}/${x}_1.fq.gz out2=${params.output}/${x}/${x}_2.fq.gz outm=${params.output}/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/${x}/${x}_phixstats.txt

   rm ${params.output}/${x}/${x}_trim*.fastq.gz
   rm ${params.output}/${x}/${x}*rmadpt.fq.gz
   
   
   """
}
