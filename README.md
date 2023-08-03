<h1 align="center">Sanibel</h1>

## What to do
The pipeline is used to analyze NGS data in fastq format from bacterial genome. It is a Nextflow version of the Flaq_amr pipeline (FL-BPHL's standard bacterial assembly pipeline with AMR detection). And some new features are added. For example, some plus analyses for specific bacteria (Legionella, Shigella, group A strep, Klebsiella, Salmonella and E.coli) and plasmid are developed. The results of these plus analyses can be found in the folder named by the bacteria name. For Neisseria and Hinfluenzae species, clonal complex and serotype prediction are automatically added to the outputs.   

## Prerequisites
Nextflow is needed. The detail of installation can be found in https://github.com/nextflow-io/nextflow.

Python3 is needed. The package "pandas" should be installed by ``` pip3 install pandas ``` if not included in your python3.

Singularity is needed. The detail of installation can be found in https://singularity-tutorial.github.io/01-installation/.

SLURM is needed.


## How to run

### Option1, your data file names directly come from Illumina output: 
1. put your data files into directory /fastqs. Your data file's name should look like "XZA22002292-XS-ASX550430-220701_S143_L001_R1_001.fastq.gz". 
2. open file "params.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output". 
3. get to the top directory of the pipeline, run 
```bash
sbatch ./sanibel_illumina.sh
```
### Option2, your data file names do not directly come from Illumina output: 
1. put your data files into directory /fastqs. Your data file's name should look like "XZA22002292_1.fastq.gz", "XZA22002292_2.fastq.gz" 
2. open file "params.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output". 
3. get into the directory of the pipeline, run 
```bash
sbatch ./sanibel.sh
```

## By Docker
By default, the pipeline uses singularity to run containers and is wrapped by SLURM. If you want to use docker to run the containers, you should use the command below:
If your data file names do not directly come from Illumina output,
```bash
sbatch ./sanibel_docker.sh
```
If your data file names directly come from Illumina output,
```bash
sbatch ./sanibel_illumina_docker.sh
```

## Version updates
    https://github.com/BPHL-Molecular/Sanibel.wiki.git
    
#### Note1: some sample data files can be found in directory /fastqs/sample_data. If you want to use these data for pipeline test, please copy them to the directory /fastqs.
#### Note2: If you want to get email notification when the pipeline running ends, please replace <EMAIL> with your email address in the line "#SBATCH --mail-user=<EMAIL>" in the batch file that you will run (namely, sanibel.sh, sanibel_illumina.sh, sanibel_docker.sh, or sanibel_illumina_docker.sh). 
#### Note3: For UF HiPerGator users, you can use "module load nextflow" to load the Nextflow in your current section, so that you donot need to install Nextflow by yourself.
#### Note4: If you want to shorten the running time, you may increase the memeory and cpus settings in the batch file that you will use. By default, memory=64gb, cpus=20, For example, you can increase the memory to 128gb by setting "#SBATCH --mem=28gb" in sanibel.sh or sanibel_illumina.sh. 