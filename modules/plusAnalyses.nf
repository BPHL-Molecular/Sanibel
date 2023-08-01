process plusAnalyses {
    input:
        val mypath
        path pyoutputs
    output:
        stdout
        //val mypath
        //path pyoutputs
        
    """
    
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    speciesid=\$(cat ${pyoutputs} | cut -d "," -f 20)
    speciesid2=\$(cat ${pyoutputs} | cut -d "," -f 1)
    echo \${speciesid}
    echo \${speciesid2}
    #the epidemiological typing of clinical and environmental isolates of Legionella pneumophila in outbreak investigations
    if [[ "\${speciesid}" == "Legionella pneumophila" || "\${speciesid2}" == "Legionella pneumophila" ]]; then
       mkdir -p ${mypath}/legsta/
       #singularity exec --cleanenv /apps/staphb-toolkit/containers/legsta_0.5.1.sif legsta ${mypath}/\${samplename}_assembly/*fasta > ${mypath}/legsta/legsta_output.txt
       singularity exec --cleanenv docker://staphb/legsta:0.5.1 legsta ${mypath}/\${samplename}_assembly/*fasta > ${mypath}/legsta/legsta_output.txt
    fi
    # determine Shigella serotype
    if [[ "\${speciesid}" == "Shigella" || "\${speciesid2}" == "Shigella" ]]; then
       mkdir -p ${mypath}/shigella/
       #singularity exec --cleanenv /apps/staphb-toolkit/containers/shigatyper_2.0.1.sif shigatyper --R1 ${params.input}/\${samplename}_1.fastq.gz --R2 ${params.input}/\${samplename}_2.fastq.gz > ${mypath}/shigella/shigella_output.txt
       singularity exec --cleanenv docker://staphb/shigatyper:2.0.1 shigatyper --R1 ${params.input}/\${samplename}_1.fastq.gz --R2 ${params.input}/\${samplename}_2.fastq.gz > ${mypath}/shigella/shigella_output.txt
    fi
    # emm type and subtype of group A strep
    if [[ "\${speciesid}" == "Streptococcus pyogenes" || "\${speciesid}" == "Streptococcus dysgalactiae" ]]; then
       mkdir -p ${mypath}/groupAstrep/
       cp ${params.input}/\${samplename}_1.fastq.gz ${mypath}/groupAstrep/
       cp ${params.input}/\${samplename}_2.fastq.gz ${mypath}/groupAstrep/
       #singularity exec -e -C -B ${mypath}/groupAstrep:/data -B ${mypath}/groupAstrep:/EMBOSS-6.6.0/emboss/.libs --pwd /data /apps/staphb-toolkit/containers/emm-typing-tool_0.0.1.sif emm_typing.py --fastq_1 ./\${samplename}_1.fastq.gz --fastq_2 ./\${samplename}_2.fastq.gz -m /db/ -o groupAstrep_output > ${mypath}/groupAstrep/groupAstrep_result.txt
       singularity exec -e -C -B ${mypath}/groupAstrep:/data -B ${mypath}/groupAstrep:/EMBOSS-6.6.0/emboss/.libs --pwd /data docker://staphb/emm-typing-tool:0.0.1 emm_typing.py --fastq_1 ./\${samplename}_1.fastq.gz --fastq_2 ./\${samplename}_2.fastq.gz -m /db/ -o groupAstrep_output > ${mypath}/groupAstrep/groupAstrep_result.txt
       rm ${mypath}/groupAstrep/\${samplename}_1.fastq.gz ${mypath}/groupAstrep/\${samplename}_2.fastq.gz
    fi
    # screen genome assemblies of Klebsiella pneumoniae and the Klebsiella pneumoniae species complex (KpSC) 
    if [[ "\${speciesid}" == "Klebsiella" || "\${speciesid2}" == "Klebsiella" ]]; then
       mkdir -p ${mypath}/klebsiella/
       cp ${mypath}/\${samplename}_assembly/\${samplename}.fasta ${mypath}/klebsiella/
       #singularity exec -e -C -B ${mypath}/klebsiella:/data --pwd /data /apps/staphb-toolkit/containers/kleborate_2.2.0.sif kleborate --resistance --kaptive --all --assemblies \${samplename}.fasta --outfile kleborate-test-out.tsv
       singularity exec -e -C -B ${mypath}/klebsiella:/data --pwd /data docker://staphb/kleborate:2.2.0 kleborate --resistance --kaptive --all --assemblies \${samplename}.fasta --outfile kleborate-test-out.tsv
       rm ${mypath}/klebsiella/\${samplename}.fasta
    fi
    # Salmonella serotype
    if [[ "\${speciesid}" == "Salmonella" || "\${speciesid2}" == "Salmonella" ]]; then
       mkdir -p ${mypath}/salmonella/
       cp ${params.input}/\${samplename}_1.fastq.gz ${mypath}/salmonella/
       cp ${params.input}/\${samplename}_2.fastq.gz ${mypath}/salmonella/
       #singularity exec -e -C -B ${mypath}/salmonella:/data --pwd /data /apps/staphb-toolkit/containers/seqsero2_1.2.1.sif SeqSero2_package.py -p 10 -t 2 -i ./\${samplename}_1.fastq.gz ./\${samplename}_2.fastq.gz
       singularity exec -e -C -B ${mypath}/salmonella:/data --pwd /data docker://staphb/seqsero2:1.2.1 SeqSero2_package.py -p 10 -t 2 -i ./\${samplename}_1.fastq.gz ./\${samplename}_2.fastq.gz
       rm ${mypath}/salmonella/\${samplename}_1.fastq.gz ${mypath}/salmonella/\${samplename}_2.fastq.gz
    fi
    # E.coli serotype finder
    if [[ "\${speciesid}" == "Escherichia coli" ]]; then
       mkdir -p ${mypath}/escherichia/
       cp ${params.input}/\${samplename}_1.fastq.gz ${mypath}/escherichia/
       cp ${params.input}/\${samplename}_2.fastq.gz ${mypath}/escherichia/
       #singularity exec -e -C -B ${mypath}/escherichia:/data --pwd /data /apps/staphb-toolkit/containers/serotypefinder_2.0.1.sif serotypefinder.py -i ./\${samplename}_1.fastq.gz ./\${samplename}_2.fastq.gz -o ./
       singularity exec -e -C -B ${mypath}/escherichia:/data --pwd /data docker://staphb/serotypefinder:2.0.1 serotypefinder.py -i ./\${samplename}_1.fastq.gz ./\${samplename}_2.fastq.gz -o ./
       rm ${mypath}/escherichia/\${samplename}_1.fastq.gz ${mypath}/escherichia/\${samplename}_2.fastq.gz
    fi
    ### find plasmid
    mkdir -p ${mypath}/plasmid/
    cp ${params.input}/\${samplename}_1.fastq.gz ${mypath}/plasmid/
    cp ${params.input}/\${samplename}_2.fastq.gz ${mypath}/plasmid/
    #singularity exec -e -C -B ${mypath}/plasmid:/data /apps/staphb-toolkit/containers/plasmidfinder_2.1.6.sif plasmidfinder.py -i /data/\${samplename}_1.fastq.gz /data/\${samplename}_2.fastq.gz -o /data/
    singularity exec -e -C -B ${mypath}/plasmid:/data docker://staphb/plasmidfinder:2.1.6 plasmidfinder.py -i /data/\${samplename}_1.fastq.gz /data/\${samplename}_2.fastq.gz -o /data/
    rm ${mypath}/plasmid/\${samplename}_1.fastq.gz ${mypath}/plasmid/\${samplename}_2.fastq.gz

    # if species is Haemophilus influenzae
    #if [[ "\${speciesid}" == "Haemophilus" || "\${speciesid2}" == "Haemophilus" ]]; then
    #   singularity exec -B ${params.output}/\${samplename}:/data docker://staphb/pmga:latest pmga --blastdir /pmga/blastdbs -o /data/pmga --force --species hinfluenzae /data/\${samplename}_assembly/\${samplename}.fasta
    #fi
    
    # if species is Neisseria
    #if [[ "\${speciesid}" == "Neisseria" || "\${speciesid2}" == "Neisseria" ]]; then
     #  singularity exec -B ${params.output}/\${samplename}:/data docker://staphb/pmga:latest pmga --blastdir /pmga/blastdbs -o /data/pmga --force --species neisseria /data/\${samplename}_assembly/\${samplename}.fasta
    #fi
    """
}
