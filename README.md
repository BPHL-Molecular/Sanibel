# Sanibel

## 1. What to do
The pipeline is a Nextflow version of the Flaq_amr pipeline. Also some plus analyses for specific bacteria (Legionella, Shigella, group A strep, Klebsiella, Salmonella and E.coli) and plasmid are developed. The results of these plus analyses can be found in the folder named by the bacteria name.  

## 2.How to run

### Option1, your data file names directly come from Illumina output: 
a. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292-FL-NDX550430-220701_S143_L001_R1_001.fastq.gz". 
b. open file "parames.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output". 
c. get to the top directory of the pipeline, run "sbatch ./sbatch_flaq_amr_nf_illumina.sh"

### Option2, your data file names do not directly come from Illumina output: 
a. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz" 
b. open file "parames.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output". 
c. get into the directory of the pipeline, run "sbatch ./sbatch_flaq_amr_nf.sh"

#### Note: some sample data files can be found in directory /fastqs/sample_data
