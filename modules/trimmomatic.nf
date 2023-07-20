process trimmomatic {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${x}"
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     

   trimmomatic PE -phred33 -trimlog ${params.output}/${x}/${x}.log ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz ${params.output}/${x}/${x}_trim_1.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/${x}/${x}_trim_2.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_2.fastq.gz ILLUMINACLIP:/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:71 TRAILING:20
   
   rm ${params.output}/${x}/${x}_unpaired_trim_*.fastq.gz
   rm ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz

   """
}
