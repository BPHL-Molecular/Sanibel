#!/usr/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=nextseq_fastq_nextflow
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20                    #This parameter shoulbe be equal to the number of samples if you want fastest running speed. However, the setting number should be less than the max cpu limit(150). 
#SBATCH --mem=20gb
#SBATCH --time=48:00:00
#SBATCH --output=flaq_amr_nf_illumina.%j.out
#SBATCH --error=flaq_amr_nf_illumina.%j.err


for f in ./fastqs/*R1*.gz
do
  cp "$f" "$(echo $f | cut -d - -f 1)_1.fastq.gz"
done
for f in ./fastqs/*R2*.gz
do
  cp "$f" "$(echo $f | cut -d - -f 1)_2.fastq.gz"
done

mkdir ./fastqs/original
mv ./fastqs/*_R1_*.gz ./fastqs/original
mv ./fastqs/*_R2_*.gz ./fastqs/original

nextflow run flaq_amr_nf.nf -params-file params.yaml

sort ./output/*/report.txt | uniq | sort -r > ./output/sum_report.txt
