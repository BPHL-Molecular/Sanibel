process unicycler {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${x}"
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   
   unicycler -1 ${params.output}/${x}/${x}_1.fq.gz -2 ${params.output}/${x}/${x}_2.fq.gz -o ${params.output}/${x}/${x}_assembly --min_fasta_length 300 --keep 1 --min_kmer_frac 0.3 --max_kmer_frac 0.9 --verbosity 2
   mv ${params.output}/${x}/${x}_assembly/assembly.fasta ${params.output}/${x}/${x}_assembly/${x}.fasta

     
   """
}
