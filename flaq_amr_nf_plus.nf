#!/usr/bin/env nextflow

/*
Note:
Before running the script, please set the parameters in the config file params.yaml
*/

//Step1:input data files
nextflow.enable.dsl=2
def L001R1Lst = []
def sampleNames = []
myDir = file("$params.input")

myDir.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}
L001R1Lst.sort()
L001R1Lst.each{
   def x = it.minus("_1.fastq.gz")
     //println x
   sampleNames.add(x)
}
//println L001R1Lst
//println sampleNames


//Step2: process the inputed data
A = Channel.fromList(sampleNames)
//A.view()

process assemble {
   input:
      val x
   output:
      //path 'xfile.txt', emit: aLook
      val "${params.output}/${x}", emit: outputpath1
      //path "${params.output}/${x}_trim_2.fastq", emit: trimR2
      
   """  
   #echo ${params.input}/${x}_1.fastq.gz >> xfile.txt
   
   mkdir -p ${params.output}/${x}
   cp ${params.input}/${x}_*.fastq.gz ${params.output}/${x}

   singularity exec --cleanenv /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz
   mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_original_fastqc.html
   mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_original_fastqc.zip
   mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_original_fastqc.html
   mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_original_fastqc.zip
   
   singularity exec --cleanenv /apps/staphb-toolkit/containers/trimmomatic_0.39.sif trimmomatic PE -phred33 -trimlog ${params.output}/${x}/${x}.log ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz ${params.output}/${x}/${x}_trim_1.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_1.fastq.gz ${params.output}/${x}/${x}_trim_2.fastq.gz ${params.output}/${x}/${x}_unpaired_trim_2.fastq.gz ILLUMINACLIP:/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:71 TRAILING:20
   rm ${params.output}/${x}/${x}_unpaired_trim_*.fastq.gz
   rm ${params.output}/${x}/${x}_1.fastq.gz ${params.output}/${x}/${x}_2.fastq.gz
   
   singularity exec --cleanenv /apps/staphb-toolkit/containers/bbtools_38.76.sif bbduk.sh in1=${params.output}/${x}/${x}_trim_1.fastq.gz in2=${params.output}/${x}/${x}_trim_2.fastq.gz out1=${params.output}/${x}/${x}_1.rmadpt.fq.gz out2=${params.output}/${x}/${x}_2.rmadpt.fq.gz ref=/bbmap/resources/adapters.fa stats=${params.output}/${x}/${x}.adapters.stats.txt ktrim=r k=23 mink=11 hdist=1 tpe tbo
   singularity exec --cleanenv /apps/staphb-toolkit/containers/bbtools_38.76.sif bbduk.sh in1=${params.output}/${x}/${x}_1.rmadpt.fq.gz in2=${params.output}/${x}/${x}_2.rmadpt.fq.gz out1=${params.output}/${x}/${x}_1.fq.gz out2=${params.output}/${x}/${x}_2.fq.gz outm=${params.output}/${x}/${x}_matchedphix.fq ref=/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats=${params.output}/${x}/${x}_phixstats.txt
   rm ${params.output}/${x}/${x}_trim*.fastq.gz
   rm ${params.output}/${x}/${x}*rmadpt.fq.gz
   
   singularity exec --cleanenv /apps/staphb-toolkit/containers/fastqc_0.11.9.sif fastqc ${params.output}/${x}/${x}_1.fq.gz ${params.output}/${x}/${x}_2.fq.gz
   mv ${params.output}/${x}/${x}_1_fastqc.html ${params.output}/${x}/${x}_1_clean_fastqc.html
   mv ${params.output}/${x}/${x}_1_fastqc.zip ${params.output}/${x}/${x}_1_clean_fastqc.zip
   mv ${params.output}/${x}/${x}_2_fastqc.html ${params.output}/${x}/${x}_2_clean_fastqc.html
   mv ${params.output}/${x}/${x}_2_fastqc.zip ${params.output}/${x}/${x}_2_clean_fastqc.zip
   
   singularity exec --cleanenv /apps/staphb-toolkit/containers/multiqc_1.8.sif multiqc ${params.output}/${x}/${x}_*_fastqc.zip -o ${params.output}/${x}
   
   
   mkdir -p ${params.output}/${x}/mash_output
   zcat ${params.output}/${x}/${x}_1.fq.gz ${params.output}/${x}/${x}_2.fq.gz | singularity exec --cleanenv /apps/staphb-toolkit/containers/mash_2.2.sif mash sketch -r -m 2 - -o ${params.output}/${x}/mash_output/${x}_sketch
   singularity exec --cleanenv /apps/staphb-toolkit/containers/mash_2.2.sif mash dist /db/RefSeqSketchesDefaults.msh ${params.output}/${x}/mash_output/${x}_sketch.msh > ${params.output}/${x}/mash_output/${x}_distances.tab
   sort -gk3 ${params.output}/${x}/mash_output/${x}_distances.tab | head > ${params.output}/${x}/mash_output/${x}_distances_top10.tab
   
   
   singularity exec --cleanenv /apps/staphb-toolkit/containers/unicycler_0.4.7.sif unicycler -1 ${params.output}/${x}/${x}_1.fq.gz -2 ${params.output}/${x}/${x}_2.fq.gz -o ${params.output}/${x}/${x}_assembly --min_fasta_length 300 --keep 1 --min_kmer_frac 0.3 --max_kmer_frac 0.9 --verbosity 2
   mv ${params.output}/${x}/${x}_assembly/assembly.fasta ${params.output}/${x}/${x}_assembly/${x}.fasta
   singularity exec --cleanenv /apps/staphb-toolkit/containers/quast_5.0.2.sif quast.py -o ${params.output}/${x}/${x}_assembly/quast_results/ ${params.output}/${x}/${x}_assembly/${x}.fasta
   """
}

process pyTask {
    input:
        val mypath
 
    output:
        //stdout
        val "${mypath}", emit: outputpath2
        path "pyoutputs.txt", emit: pyoutputs
 
    """
    #!/usr/bin/env python3
    import sys
    import os
    import subprocess
    import argparse
    import datetime
    import fnmatch
    import re
    import pandas as pd
    
    items = "${mypath}".strip().split("/")
    #print(items[-1])
    filepath1 = "${mypath}"+"/mash_output/"+items[-1]+"_distances_top10.tab"
    with open(filepath1, 'r') as mash:
        top_hit = mash.readline()
        top_hit = str(top_hit)
        cells = top_hit.split()
        acell = cells[0]
        subcells = re.split('-.-|_subsp',acell)
        subsubcells = subcells[1].split('_')
        genus = subsubcells[0]
        species = subsubcells[1]
        distance = cells[2]
        accession = top_hit.split("-")[5]
    #print(accession)
    
    filepath2 = "${mypath}"+"/"+items[-1]+"_assembly/quast_results/report.tsv"
    df = pd.read_table(filepath2, sep="\t")
    assem = list(df.columns)[1]
    #print(assem)
    contigs = df[assem][12].astype(int)
    long_contig = df[assem][13].astype(int)
    n50 = df[assem][16].astype(int)
    l50 = df[assem][18].astype(int)
    genome = df[assem][14].astype(int)
    gc = df[assem][15].astype(int)
    #print(genome)
    

    f = open("pyoutputs.txt", "w")
    f.write(str(genus)+","+str(species)+","+str(distance)+","+str(accession)+","+str(assem)+","+str(contigs)+","+str(long_contig)+","+str(n50)+","+str(l50)+","+str(genome)+","+str(gc))
    f.close
    
    """
}

process readssum {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    #echo ${mypath}
    genome=\$(cat ${pyoutputs} | cut -d "," -f 10)
    #echo \${genome}
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    #echo \${samplename}
    singularity exec --cleanenv /apps/staphb-toolkit/containers/lyveset_1.1.4f.sif run_assembly_shuffleReads.pl -gz ${mypath}/\${samplename}_1.fq.gz ${mypath}/\${samplename}_2.fq.gz > ${mypath}/\${samplename}_clean_shuffled.fq.gz
    singularity exec --cleanenv /apps/staphb-toolkit/containers/lyveset_1.1.4f.sif run_assembly_readMetrics.pl ${mypath}/\${samplename}_clean_shuffled.fq.gz -e \${genome} > ${mypath}/\${samplename}_readMetrics.txt
    """
}


process pyTask2 {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    #!/usr/bin/env python3
    items = "${mypath}".strip().split("/")
    filepath = "${mypath}/"+items[-1]+"_readMetrics.txt"
    with open(filepath, 'r') as metrics:
        firstline = metrics.readline()
        firstline = firstline.rstrip()
        firstline = firstline.split()
        secondline = metrics.readline()
        secondline = secondline.rstrip()
        secondline = secondline.split()
        rm = dict(zip(firstline, secondline))
        avg_read_len = rm['avgReadLength']
        avg_qual = rm['avgQuality']
        num_read = rm['numReads']
        cov = rm['coverage']
    #print(cov)
    
    #f = open("pyoutputs2.txt", "w")
    f = open("${pyoutputs}", "a")
    f.write(","+str(avg_read_len)+","+str(avg_qual)+","+str(num_read)+","+str(cov))
    f.close
    """
}

process annotate {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    samplename=\$(echo ${mypath} | rev | cut -d "/" -f 1 | rev)
    genus=\$(cat ${pyoutputs} | cut -d "," -f 1)
    species=\$(cat ${pyoutputs} | cut -d "," -f 2)
    echo \${genus}
    singularity exec --cleanenv /apps/staphb-toolkit/containers/prokka_1.14.5.sif prokka --genus \${genus} --species \${species} --strain \${samplename} --outdir ${mypath}/\${samplename}_assembly/prokka --prefix \${samplename} --force --compliant --locustag \${genus} ${mypath}/\${samplename}_assembly/\${samplename}.fasta
    singularity exec --cleanenv -B /tmp:/tmp -B ${mypath}/\${samplename}_assembly/:/data /apps/staphb-toolkit/containers/ncbi-amrfinderplus_3.10.1.sif amrfinder -n /data/\${samplename}.fasta --plus -o /data/\${samplename}_amrfinderplus_report.tsv
    
    singularity exec --bind ${mypath}/\${samplename}_assembly:/data --pwd /data --cleanenv /apps/staphb-toolkit/containers/mlst_2.19.0.sif mlst ${mypath}/\${samplename}_assembly/\${samplename}.fasta --nopath > ${mypath}/\${samplename}_assembly/\${samplename}.mlst
    
    mkdir -p ${mypath}/kraken_out/
    singularity exec --cleanenv /apps/staphb-toolkit/containers/kraken2_2.0.8-beta.sif kraken2 --db /kraken2-db/minikraken2_v1_8GB/ --use-names --report ${mypath}/kraken_out/\${samplename}.report --output ${mypath}/kraken_out/\${samplename}_kraken.out --paired ${params.input}/\${samplename}_1.fastq.gz ${params.input}/\${samplename}_2.fastq.gz
    
    mkdir -p assemblies
    mkdir -p annotations
    mkdir -p amrfinder_results
    cp ${mypath}/\${samplename}_assembly/\${samplename}.fasta assemblies/
    cp ${mypath}/\${samplename}_assembly/prokka/\${samplename}.gff annotations/
    cp ${mypath}/\${samplename}_assembly/\${samplename}_amrfinderplus_report.tsv amrfinder_results/ 
    """
}

process pyTask3 {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    #!/usr/bin/env python3
    items = "${mypath}".strip().split("/")
    filepath1="${mypath}/"+items[-1]+"_assembly/prokka/"+items[-1]+".txt"
    filepath2="${mypath}/"+items[-1]+"_assembly/"+items[-1]+".mlst"
    filepath3="${mypath}/kraken_out/"+items[-1]+".report"
    filepath4="${mypath}/report.txt"
    
    cds = ''
    with open(filepath1, 'r') as genes:
        for line in genes:
            line = line.rstrip()
            content = line.split()
            if content[0] == 'CDS:':
                cds = content[1]
    with open(filepath2,'r') as mlst:
       for line in mlst:
           out = line.rstrip().split()
           scheme = out[1]
           st = out[2]
    with open(filepath3, 'r') as kreport:
        lines = kreport.readlines()
        for l in lines:
            l_parse = l.lstrip().rstrip().split("\t")
            percent = l_parse[0]
            tax_level = l_parse[3]
            tax = l_parse[5].lstrip()
            if tax_level == 'S':
                break
    #print(cds+"\t"+st+"\t"+tax)
    
    f = open("${pyoutputs}", "a")
    f.write(","+str(cds)+","+str(scheme)+","+str(st)+","+str(percent)+","+str(tax))
    f.close
    
    #f=open("${pyoutputs}", "r")
    #print(f.read())
    
    """
}

process pyTask4 {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    $/
    #!/usr/bin/env python3
    items = "${mypath}".strip().split("/")
    #parentpath = "/".join(items[:-1])
    
    with open("${pyoutputs}", "r") as aline:
        for line in aline:
            cells=line.rstrip().split(",")
            results=[items[-1],cells[0]+'_'+cells[1],cells[3],cells[2],cells[19],cells[18],cells[16],cells[17],cells[13],cells[11],cells[12],cells[14],cells[5],cells[6],cells[7],cells[8],cells[9],cells[10],cells[15]]
            print(results)

    report = open("${mypath}"+"/report.txt", 'w')
    header = ['sampleID', 'speciesID_mash', 'nearest_neighbor_mash', 'mash_distance', 'speciesID_kraken', 'kraken_percent', 'mlst_scheme', 'mlst_st', 'num_clean_reads', 'avg_readlength', 'avg_read_qual', 'est_coverage', 'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length', 'gc_content', 'annotated_cds'] 
    report.write("\t".join(header))
    report.write('\n')
    report.write('\t'.join(results))
    report.close()
    /$
}

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
       singularity exec --cleanenv /apps/staphb-toolkit/containers/legsta_0.5.1.sif legsta ${mypath}/\${samplename}_assembly/*fasta > ${mypath}/legsta/legsta_output.txt
    fi
    # determine Shigella serotype
    if [[ "\${speciesid}" == "Shigella" || "\${speciesid2}" == "Shigella" ]]; then
       mkdir -p ${mypath}/shigella/
       singularity exec --cleanenv /apps/staphb-toolkit/containers/shigatyper_2.0.1.sif shigatyper --R1 ${params.input}/\${samplename}_1.fastq.gz --R2 ${params.input}/\${samplename}_2.fastq.gz > ${mypath}/shigella/shigella_output.txt
    fi
    # emm type and subtype of group A strep
    if [[ "\${speciesid}" == "Streptococcus pyogenes" || "\${speciesid}" == "Streptococcus dysgalactiae" ]]; then
       mkdir -p ${mypath}/groupAstrep/
       cp ${params.input}/\${samplename}_1.fastq.gz ${mypath}/groupAstrep/
       cp ${params.input}/\${samplename}_2.fastq.gz ${mypath}/groupAstrep/
       singularity exec -e -C -B ${mypath}/groupAstrep:/data -B ${mypath}/groupAstrep:/EMBOSS-6.6.0/emboss/.libs --pwd /data /apps/staphb-toolkit/containers/emm-typing-tool_0.0.1.sif emm_typing.py --fastq_1 ./\${samplename}_1.fastq.gz --fastq_2 ./\${samplename}_2.fastq.gz -m /db/ -o groupAstrep_output > ${mypath}/groupAstrep/groupAstrep_result.txt
       rm ${mypath}/groupAstrep/\${samplename}_1.fastq.gz ${mypath}/groupAstrep/\${samplename}_2.fastq.gz
    fi
    # screen genome assemblies of Klebsiella pneumoniae and the Klebsiella pneumoniae species complex (KpSC) 
    if [[ "\${speciesid}" == "Klebsiella" || "\${speciesid2}" == "Klebsiella" ]]; then
       mkdir -p ${mypath}/klebsiella/
       cp ${mypath}/\${samplename}_assembly/\${samplename}.fasta ${mypath}/klebsiella/
       singularity exec -e -C -B ${mypath}/klebsiella:/data --pwd /data /apps/staphb-toolkit/containers/kleborate_2.2.0.sif kleborate --resistance --kaptive --all --assemblies \${samplename}.fasta --outfile kleborate-test-out.tsv
       rm ${mypath}/klebsiella/\${samplename}.fasta
    fi
    # Salmonella serotype
    if [[ "\${speciesid}" == "Salmonella" || "\${speciesid2}" == "Salmonella" ]]; then
       mkdir -p ${mypath}/salmonella/
       cp ${params.input}/\${samplename}_1.fastq.gz ${mypath}/salmonella/
       cp ${params.input}/\${samplename}_2.fastq.gz ${mypath}/salmonella/
       singularity exec -e -C -B ${mypath}/salmonella:/data --pwd /data /apps/staphb-toolkit/containers/seqsero2_1.2.1.sif SeqSero2_package.py -p 10 -t 2 -i ./\${samplename}_1.fastq.gz ./\${samplename}_2.fastq.gz
       rm ${mypath}/salmonella/\${samplename}_1.fastq.gz ${mypath}/salmonella/\${samplename}_2.fastq.gz
    fi
    # E.coli serotype finder
    if [[ "\${speciesid}" == "Escherichia coli" ]]; then
       mkdir -p ${mypath}/escherichia/
       cp ${params.input}/\${samplename}_1.fastq.gz ${mypath}/escherichia/
       cp ${params.input}/\${samplename}_2.fastq.gz ${mypath}/escherichia/
       singularity exec -e -C -B ${mypath}/escherichia:/data --pwd /data /apps/staphb-toolkit/containers/serotypefinder_2.0.1.sif serotypefinder.py -i ./\${samplename}_1.fastq.gz ./\${samplename}_2.fastq.gz -o ./
       rm ${mypath}/escherichia/\${samplename}_1.fastq.gz ${mypath}/escherichia/\${samplename}_2.fastq.gz
    fi
    ### find plasmid
    mkdir -p ${mypath}/plasmid/
    cp ${params.input}/\${samplename}_1.fastq.gz ${mypath}/plasmid/
    cp ${params.input}/\${samplename}_2.fastq.gz ${mypath}/plasmid/
    singularity exec -e -C -B ${mypath}/plasmid:/data /apps/staphb-toolkit/containers/plasmidfinder_2.1.6.sif plasmidfinder.py -i /data/\${samplename}_1.fastq.gz /data/\${samplename}_2.fastq.gz -o /data/
    rm ${mypath}/plasmid/\${samplename}_1.fastq.gz ${mypath}/plasmid/\${samplename}_2.fastq.gz
    """
}

workflow {
    assemble(A) | pyTask | readssum | pyTask2 |  annotate| pyTask3 | pyTask4 | plusAnalyses | view
}
