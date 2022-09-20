<h1 align="center">Sanibel</h1>

## What to do
The pipeline is a Nextflow version of the Flaq_amr pipeline (FL-BPHL's standard bacterial assembly pipeline with AMR detection). Also some plus analyses for specific bacteria (Legionella, Shigella, group A strep, Klebsiella, Salmonella and E.coli) and plasmid are developed. The results of these plus analyses can be found in the folder named by the bacteria name.  

## Prerequisites
Nextflow is needed. The detail of installation can be found in https://github.com/nextflow-io/nextflow.
Python3 is needed.
Singularity is also needed. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/.
In addition, the below docker container images are needed in the pipeline. These images should be downloaded to the directory /apps/staphb-toolkit/containers/ in your local computer. You can find them from StaPH-B/docker-builds (https://github.com/StaPH-B/docker-builds).
1. trimmomatic_0.39.sif
2. bbtools_38.76.sif
3. fastqc_0.11.9.sif
4. multiqc_1.8.sif
5. mash_2.2.sif
6. unicycler_0.4.7.sif
7. quast_5.0.2.sif
8. lyveset_1.1.4f.sif
9. prokka_1.14.5.sif
10. ncbi-amrfinderplus_3.10.1.sif
11. mlst_2.19.0.sif
12. kraken2_2.0.8-beta.sif
13. legsta_0.5.1.sif
14. shigatyper_2.0.1.sif
15. emm-typing-tool_0.0.1.sif
16. kleborate_2.2.0.sif
17. seqsero2_1.2.1.sif
18. serotypefinder_2.0.1.sif
19. plasmidfinder_2.1.6.sif
## How to run

### Option1, your data file names directly come from Illumina output: 
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292-FL-NDX550430-220701_S143_L001_R1_001.fastq.gz". 
2. open file "parames.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output". 
3. get to the top directory of the pipeline, run "sbatch ./sbatch_flaq_amr_nf_illumina.sh"

### Option2, your data file names do not directly come from Illumina output: 
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz" 
2. open file "parames.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output". 
3. get into the directory of the pipeline, run "sbatch ./sbatch_flaq_amr_nf.sh"

#### Note: some sample data files can be found in directory /fastqs/sample_data
