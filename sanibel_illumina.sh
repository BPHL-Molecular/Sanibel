#!/usr/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=sanibel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20                    #This parameter shoulbe be equal to the number of samples if you want fastest running speed. However, the setting number should be less than the max cpu limit(150). 
#SBATCH --mem=200gb
#SBATCH --time=48:00:00
#SBATCH --output=sanibel.%j.out
#SBATCH --error=sanibel.%j.err
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=FAIL,END

module load nextflow
APPTAINER_CACHEDIR=./
export APPTAINER_CACHEDIR


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

singularity exec docker://staphb/mlst:2.23.0 cp /mlst-2.23.0/db/pubmlst/neisseria/neisseria.txt ./
singularity exec  docker://staphb/mlst:2.23.0 cp /mlst-2.23.0/db/pubmlst/hinfluenzae/hinfluenzae.txt ./

nextflow run flaq_amr_plus2.nf -params-file params.yaml

sort ./output/*/report.txt | uniq > ./output/sum_report.txt
sed -i '/sampleID\tspeciesID/d' ./output/sum_report.txt
sed -i '1i sampleID\tspeciesID_mash\tnearest_neighb_mash\tmash_distance\tspeciesID_kraken\tkraken_percent\tmlst_scheme\tmlst_st\tmlst_cc\tpmga_species\tserotype\tnum_clean_reads\tavg_readlength\tavg_read_qual\test_coverage\tnum_contigs\tlongest_contig\tN50\tL50\ttotal_length\tgc_content\tannotated_cds' ./output/sum_report.txt
rm ./neisseria.txt
rm ./hinfluenzae.txt

rm ./fastqs/*fastq.gz
mv ./fastqs/original/*fastq.gz ./fastqs/
rm -r ./fastqs/original

mv ./*.out ./output
mv ./*err ./output

dt=$(date "+%Y%m%d%H%M%S")
mv ./output ./output-$dt
#mv ./work ./work-$dt
rm -r ./work
rm -r ./cache

