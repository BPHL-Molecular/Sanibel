process assemble {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${params.output}/${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """     
   mkdir -p ${params.output}/${x}
   cp ${params.input}/${x}_*.fastq.gz ${params.output}/${x}

   ##### pull container from docker hub, and if multiple containers in a process
   singularity exec --cleanenv docker://staphb/fastqc:0.11.9 fastqc ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz

   mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_original_fastqc.html
   mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_original_fastqc.zip
   mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_original_fastqc.html
   mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_original_fastqc.zip
   
   singularity exec --cleanenv docker://staphb/trimmomatic:0.39 trimmomatic PE -phred33 -trimlog ${params.output}/${x}/${x}.log ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz ${params.output}/${x}/${x}_trim_1.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/${x}/${x}_trim_2.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_2.fastq.gz ILLUMINACLIP:/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:71 TRAILING:20
   
   rm ${params.output}/${x}/${x}_unpaired_trim_*.fastq.gz
   rm ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz

   singularity exec --cleanenv docker://staphb/bbtools:38.76 bbduk.sh in1=${params.output}/${x}/${x}_trim_1.fastq.gz in2=${params.output}/${x}/${x}_trim_2.fastq.gz out1=${params.output}/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo

   singularity exec --cleanenv docker://staphb/bbtools:38.76 bbduk.sh in1=${params.output}/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/${x}/${x}_1.fq.gz out2=${params.output}/${x}/${x}_2.fq.gz outm=${params.output}/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/${x}/${x}_phixstats.txt

   rm ${params.output}/${x}/${x}_trim*.fastq.gz
   rm ${params.output}/${x}/${x}*rmadpt.fq.gz
   
   singularity exec --cleanenv docker://staphb/fastqc:0.11.9 fastqc ${params.output}/${x}/${x}_1.fq.gz ${params.output}/${x}/${x}_2.fq.gz

   mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_clean_fastqc.html
   mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_clean_fastqc.zip
   mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_clean_fastqc.html
   mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_clean_fastqc.zip
   
   singularity exec --cleanenv docker://staphb/multiqc:1.8 multiqc ${params.output}/${x}/${x}_*_fastqc.zip -o ${params.output}/${x}
     
   mkdir -p ${params.output}/${x}/mash_output
   zcat ${params.output}/${x}/${x}_1.fq.gz ${params.output}/${x}/${x}_2.fq.gz | singularity exec --cleanenv /apps/staphb-toolkit/containers/mash_2.2.sif mash sketch -r -m 2 - -o ${params.output}/${x}/mash_output/${x}_sketch
   
   singularity exec --cleanenv docker://staphb/mash:2.2 mash dist /db/RefSeqSketchesDefaults.msh ${params.output}/${x}/mash_output/${x}_sketch.msh > ${params.output}/${x}/mash_output/${x}_distances.tab
   sort -gk3 ${params.output}/${x}/mash_output/${x}_distances.tab | head > ${params.output}/${x}/mash_output/${x}_distances_top10.tab
     
   singularity exec --cleanenv docker://staphb/unicycler:0.4.7 unicycler -1 ${params.output}/${x}/${x}_1.fq.gz -2 ${params.output}/${x}/${x}_2.fq.gz -o ${params.output}/${x}/${x}_assembly --min_fasta_length 300 --keep 1 --min_kmer_frac 0.3 --max_kmer_frac 0.9 --verbosity 2
   mv ${params.output}/${x}/${x}_assembly/assembly.fasta ${params.output}/${x}/${x}_assembly/${x}.fasta
   if [[ ${params.isNeisseria} =~ yes ]];then
      singularity exec -B ${params.output}/${x}:/data docker://staphb/pmga:latest pmga --blastdir /pmga/blastdbs -o /data/pmga --force --species neisseria /data/${x}_assembly/${x}.fasta
   fi
   if [[ ${params.isHinfluenzae} =~ yes ]];then
      singularity exec -B ${params.output}/${x}:/data docker://staphb/pmga:latest pmga --blastdir /pmga/blastdbs -o /data/pmga --force --species hinfluenzae /data/${x}_assembly/${x}.fasta
   fi
   singularity exec --cleanenv docker://staphb/quast:5.0.2 quast.py -o ${params.output}/${x}/${x}_assembly/quast_results/ ${params.output}/${x}/${x}_assembly/${x}.fasta



   """
}
