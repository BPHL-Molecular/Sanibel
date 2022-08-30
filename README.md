# Sanibel
Nextflow version of the Flaq_amr pipeline

1.How to run
1) If you data file names directly come from Illumina output:
   a. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292-FL-NDX550430-220701_S143_L001_R1_001.fastq.gz"
   b. open file "parames.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output".
   c. get to the top directory of the pipeline, run "sbatch ./sbatch_flaq_amr_nf_illumina.sh"

2) If you data file names do not directly come from Illumina output:
   a. put your data files into directory /fastqs. Your data file's name should look like "JBS22002292_1.fastq.gz", "JBS22002292_2.fastq.gz"
   b. open file "parames.yaml", set the two parameters absolute paths. They should be ".../.../fastqs" and ".../.../output".
   c. get into the directory of the pipeline, run "sbatch ./sbatch_flaq_amr_nf.sh"
   
Note: some sample data files can be found in directory /sample_data
