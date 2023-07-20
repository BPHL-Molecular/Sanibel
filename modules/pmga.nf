process pmga {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${x}"
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   
   if [[ ${params.isNeisseria} =~ yes ]];then
      #singularity exec -B ${params.output}/${x}:/data docker://staphb/pmga:latest pmga --blastdir /pmga/blastdbs -o /data/pmga --force --species neisseria /data/${x}_assembly/${x}.fasta
      pmga --blastdir /pmga/blastdbs -o ${params.output}/${x}/pmga --force --species neisseria ${params.output}/${x}/${x}_assembly/${x}.fasta
   fi
   if [[ ${params.isHinfluenzae} =~ yes ]];then
      #singularity exec -B ${params.output}/${x}:/data docker://staphb/pmga:latest pmga --blastdir /pmga/blastdbs -o /data/pmga --force --species hinfluenzae /data/${x}_assembly/${x}.fasta
      pmga --blastdir /pmga/blastdbs -o ${params.output}/${x}/pmga --force --species hinfluenzae ${params.output}/${x}/${x}_assembly/${x}.fasta
   fi
     
   """
}
