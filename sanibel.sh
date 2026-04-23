#!/usr/bin/bash
#SBATCH --account=bphl-umbrella
#SBATCH --qos=bphl-umbrella
#SBATCH --job-name=sanibel
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200gb
#SBATCH --time=48:00:00
#SBATCH --output=sanibel.%j.out
#SBATCH --error=sanibel.%j.err
#SBATCH --mail-user=<EMAIL>
#SBATCH --mail-type=FAIL,END

module load conda nextflow apptainer
conda activate SANIBEL

# Path to container image cache directory
export NXF_APPTAINER_CACHEDIR=/path/to/singularity/cache

# Run pipeline
nextflow run sanibel.nf -profile apptainer -params-file params.yaml

# Rename output directory with timestamp
dt=$(date "+%Y%m%d%H%M%S")
output_dir=$(grep '^output:' params.yaml | sed 's/output:[[:space:]]*//' | tr -d '"')
mv "$output_dir" "${output_dir}-${dt}"

# Cleanup (disabled for troubleshooting runs)
#rm -rf ./work ./cache
