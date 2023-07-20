process quast {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${params.output}/${x}"
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   
   quast.py -o ${params.output}/${x}/${x}_assembly/quast_results/ ${params.output}/${x}/${x}_assembly/${x}.fasta

     
   """
}
