<p align="center">
  <img src="assets/sanibel_pipeline_logo_v2.svg" alt="Sanibel logo" width="600"/>
</p>

<p align="center">
  <em>⚠️ For research use only. Results were obtained by procedures that were not CLIA validated.</em>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Pipeline-Sanibel%202.1-blue?style=plastic" />
  <img src="https://img.shields.io/badge/Nextflow-≥23.04-brightgreen?style=plastic&logo=nextflow" />
  <img src="https://img.shields.io/badge/Python-3.10+-yellow?style=plastic&logo=python" />
  <img src="https://img.shields.io/badge/License-Apache%202.0-red?style=plastic" />
</p>

## 🦠🧬 Overview

Sanibel is Florida BPHL's Nextflow bacterial whole-genome sequencing (WGS) analysis pipeline. It performs quality control, *de novo* assembly, taxonomic classification, species identification, sequence typing and antimicrobial resistance (AMR) detection on paired-end Illumina short reads.

Species identification uses a candidate-pool design: Mash, Kraken2 and 16S rRNA BLAST nominate a ranked pool of candidate species, multiple RefSeq genomes are downloaded per candidate and skani confirms the species by whole-genome ANI. skani is the arbiter for the reported organism and the contamination flag. A confident call requires ANI ≥ 95% and alignment fraction ≥ 50% for a positive ID.

Species-specific typing modules run automatically based on species identification results: *Legionella pneumophila* (Legsta), *Klebsiella* (Kleborate), *Shigella* (ShigaTyper), *Streptococcus pyogenes/dysgalactiae* (EMM typing), *Salmonella* (SeqSero2), *E. coli* (SerotypeFinder), *Streptococcus pneumoniae* (SeroBA), *Pseudomonas aeruginosa* (pasty), *Acinetobacter baumannii* (Kaptive), *Vibrio parahaemolyticus* (Kaptive), *Listeria monocytogenes* (LisSero) and *Neisseria meningitidis*/*Haemophilus influenzae* (PMGA). BMGAP2 provides enhanced AMR and antigen analysis for *Neisseria meningitidis* and *Haemophilus influenzae*. PlasmidFinder runs on all samples.

### ⚙️ Dependencies

- **Nextflow** 23.04–26.x - [installation guide](https://github.com/nextflow-io/nextflow)
- **Apptainer/Singularity** - [installation guide](https://apptainer.org/docs/user/latest/)
- **Conda** - [installation guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- **SLURM** workload manager (required for HiPerGator; otherwise not required)

All bioinformatics tools run inside containers, no additional software installation is required.

### 💻 Resource Requirements

Sanibel can run on any system with Nextflow and Apptainer installed, but is **strongly recommended to run on an HPC environment**. The most resource intensive step is **Unicycler** (*de novo* assembly via SPAdes), which runs once per sample.

- **CPUs:** 20 recommended (runs 2 assemblies in parallel); minimum 4
- **RAM:** 64 GB recommended; minimum 16 GB
- **Disk:** ~10 GB per sample (input + output); ~8.5 GB for the Kraken2 database (if using minikraken or standard-8 Kraken2 databases)

> **Running locally with fewer CPUs:** Nextflow will still run but your machine will be oversubscribed during assembly. Each Unicycler job requests 10 CPUs by default, on a machine with fewer cores the OS will time-share threads and assembly will complete more slowly but will not fail. You can lower the `cpus` value for `unicycler` in `nextflow.config` to match your hardware.

**Estimated runtime** (12 samples, 20 CPUs, HPC): ~2 hours total.

### 🛠️ Setup

#### 1. Clone this repository and enter the repository directory

```bash
$ git clone https://github.com/BPHL-Molecular/Sanibel
$ cd Sanibel/
```

#### 2. Create the conda environment

```bash
$ conda env create -f environment.yaml
$ conda activate SANIBEL
```

#### 3. Configure params.yaml

Edit `params.yaml` and set paths for your environment:

```yaml
# Input / Output absolute path, no trailing slash "/"
input:  "/full/path/to/fastqs"
output: "/full/path/to/output"

# non-Florida-BPHL users: uncomment and set your BMGAP2 analysis_scripts directory
#bmgap2_db:  "/full/path/to/bmgap2/analysis_scripts"

# non-Florida-BPHL users: uncomment and set your Kraken2 database directory
#kraken_db:  "/full/path/to/kraken2/database"
```

> **Florida BPHL users:** only `input` and `output` need to be set. All other paths are pre-configured.
> **Non-Florida BPHL users:** uncomment `bmgap2_db` and `kraken_db` by removing the leading `#` and set their paths.

#### 4. Configure sanibel.sh

> At Florida BPHL we use **Apptainer** on HiPerGator for containerization. `sanibel.sh` is pre-configured for SLURM + Apptainer and is the recommended submission method for HiPerGator users.

Add your email address for job notifications and set `NXF_APPTAINER_CACHEDIR` to your image cache directory:

```bash
#SBATCH --mail-user=your@email.gov
export NXF_APPTAINER_CACHEDIR=/path/to/apptainer/cache

```

#### 5. Kraken2 and BMGAP2 Setup for Non-Florida-BPHL Users

> **Florida BPHL users:** Kraken2 and BMGAP2 are already installed and configured on the cluster. The `kraken_db` and `bmgap2_db` paths in `params.yaml` are pre-set. Skip this section entirely.

For non-Florida BPHL users, a Kraken2 database must be downloaded before running the pipeline. These are available [here](https://benlangmead.github.io/aws-indexes/k2). Once downloaded, set the location path in the `params.yaml` file.

[BMGAP2](https://github.com/CDCgov/BMGAP2) (*Bacterial Meningitis Genome Analysis Pipeline 2*) runs automatically on *Neisseria meningitidis* and *Haemophilus influenzae* samples. The custom python scripts check the MLST scheme internally and skip any sample that is not *N. meningitidis* or *H. influenzae*.

For non-Florida BPHL users, BMGAP2 must be installed before running the pipeline. Its scripts run directly on the host (not inside a container) and are invoked by the three `bmgap2_*` Nextflow modules. Follow the [BMGAP2 installation instructions](https://github.com/CDCgov/BMGAP2) to clone the repository and build all required databases, then set `bmgap2_db` in `params.yaml` to the `analysis_scripts` directory.

### How to Run

Place input FASTQ files in the directory specified by `params.input`. Both Illumina native (`SAMPLE_S1_L001_R1_001.fastq.gz`) and simplified (`SAMPLE_1.fastq.gz`) naming conventions are supported.

### 🐊 HiPerGator Usage
```bash
sbatch sanibel.sh
```

### ⚡ Local Usage
```bash
nextflow run sanibel.nf -profile apptainer -params-file params.yaml
```

### Workflow Diagram

```mermaid
flowchart LR
    IN[Paired FASTQ] --> QC["Quality control<br/>FastQC · Trimmomatic · BBTools"]
    QC --> ASM["Assembly<br/>Unicycler · Quast"]
    ASM --> SID["Species ID<br/>Mash · Kraken2 · 16S rRNA BLAST<br/>confirmed by skani ANI"]
    SID --> TYP["Annotation and Typing<br/>Prokka · MLST · Serotyping"]

    ASM --> AMR["AMR and plasmids<br/>AMRFinder · PlasmidFinder"]
    QC --> AMR

    QC --> REP[summary_report]
    ASM --> REP
    SID --> REP
    TYP --> REP
    AMR --> REP

    REP --> OUT["sum_report.txt<br/>nm_sum_report.txt<br/>hi_sum_report.txt"]
    REP --> MQ["multiqc_global<br/>sanibel_report.html"]

    style SID fill:#9f9,stroke:#333,color:#000
    style REP fill:#f96,stroke:#333,stroke-width:2px,color:#000
    style OUT fill:#f96,stroke:#333,color:#000
    style MQ fill:#f96,stroke:#333,color:#000
```

For how the species ID vote, candidate pool and skani ANI confirmation fit together, see [docs/species-id.md](docs/species-id.md).

### 🧩 Modules

Sanibel is made possible thanks to the following tools:

<small>

**Quality Control** - [FastQC](https://github.com/s-andrews/FastQC) · [Trimmomatic](https://github.com/usadellab/Trimmomatic) · [BBTools](https://github.com/bbushnell/BBTools) · [MultiQC](https://github.com/MultiQC/MultiQC) · [Lyveset](https://github.com/lskatz/lyve-SET)

**Assembly & Annotation** - [Unicycler](https://github.com/rrwick/Unicycler) · [QUAST](https://github.com/ablab/quast) · [Prokka](https://github.com/tseemann/prokka)

**Typing & Classification** - [Mash](https://github.com/marbl/Mash) · [Kraken2](https://github.com/DerrickWood/kraken2) · [BLAST](https://blast.ncbi.nlm.nih.gov) · [NCBI Datasets](https://github.com/ncbi/datasets) · [skani](https://github.com/bluenote-1577/skani) · [MLST](https://github.com/tseemann/mlst)

**Species-Specific Serotyping** - [Legsta](https://github.com/MDU-PHL/legsta) · [Kleborate](https://github.com/klebgenomics/Kleborate) · [ShigaTyper](https://github.com/CFSAN-Biostatistics/shigatyper) · [emm-typing-tool](https://github.com/ukhsa-collaboration/emm-typing-tool) · [SeqSero2](https://github.com/denglab/SeqSero2) · [SerotypeFinder](https://bitbucket.org/genomicepidemiology/serotypefinder) · [PMGA](https://github.com/CDCgov/PMGA) · [BMGAP2](https://github.com/CDCgov/BMGAP2) · [SeroBA](https://github.com/sanger-pathogens/seroba) · [pasty](https://github.com/rpetit3/pasty) · [Kaptive](https://github.com/klebgenomics/Kaptive) · [LisSero](https://github.com/MDU-PHL/LisSero)

**AMR & Mobile Genetic Elements** - [AMRFinderPlus](https://github.com/ncbi/amr) · [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder)

</small>

### 📁 Output

All results are written to `params.output/<sample_id>/`. Depending on which species are in the run, up to three summary files are written to `params.output/`:

| File | Samples | Cols | Key fields |
|------|---------|------|------------|
| `sum_report.txt` | All | 31 | ID · species (skani ANI, Mash, Kraken) · 16S top hit · skani ANI/reference · species-ID QC · contamination flag · MLST scheme/ST · serotype · QC metrics (reads, coverage, assembly stats, GC, CDS) · assembly QC · AMR gene symbols/subclasses |
| `nm_sum_report.txt` | *N. meningitidis* only | 26 | ID · PMGA serogroup · BMGAP2 AMR alleles/phenotypes · vaccine antigen coverage (4CMenB) |
| `hi_sum_report.txt` | *H. influenzae* only | 22 | ID · PMGA capsule type · BMGAP2 AMR alleles/phenotypes |


### 🤝 Contributing
We welcome contributions to make Sanibel better! Feel free to open issues or submit pull requests to suggest any additional features or enhancements!

### 📧 Contact
**Email**: bphl-sebioinformatics@flhealth.gov

### ⚖️ License
Sanibel is licensed under the [Apache License](https://github.com/BPHL-Molecular/Sanibel/blob/main/LICENSE).

[BMGAP2](https://github.com/CDCgov/BMGAP2) is developed by the CDC and is also distributed under the Apache License 2.0. A copy is included in [`licenses/BMGAP2_LICENSE`](licenses/BMGAP2_LICENSE).
