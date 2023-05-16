<h1 align="center">Sanibel</h1>

## What to do
The pipeline is a Nextflow version of the Flaq_amr pipeline (FL-BPHL's standard bacterial assembly pipeline with AMR detection). Also some plus analyses for specific bacteria (Legionella, Shigella, group A strep, Klebsiella, Salmonella and E.coli) and plasmid are developed. The results of these plus analyses can be found in the folder named by the bacteria name. In addition, for Neisseria and Hinfluenzae, clonal complex and serotype prediction are added to the analysis.   

## Prerequisites
Nextflow is needed. The detail of installation can be found in https://github.com/nextflow-io/nextflow.

Python3 is needed.

Singularity is also needed. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/.


## How to run

### Option1, your data file names directly come from Illumina output: 
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292-FL-NDX550430-220701_S143_L001_R1_001.fastq.gz". 
2. open file "parames.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output". 
3. get to the top directory of the pipeline, run "sbatch ./sbatch_flaq_amr_nf_illumina.sh"

### Option2, your data file names do not directly come from Illumina output: 
1. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz" 
2. open file "parames.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output". 
3. get into the directory of the pipeline, run "sbatch ./sbatch_flaq_amr_nf.sh"

## For Neisseria or Hinfluenzae analysis

For Neisseria analysis, the value of the parameter "isNeisseria" is set as "yes"

For Hinfluenzae analysis, the value of the parameter "isHinfluenzae" is set as "yes"

By default the values of the above two parameters are both "no". Although Neisseria or Hinfluenzae can be analyzed using their default settings (namely "no"), clonal complex and serotype prediction features will not work.   


#### Note: some sample data files can be found in directory /fastqs/sample_data
