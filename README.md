<h1 align="center">Sanibel - Bacterial WGS Analysis Pipeline</h1>

<p align="center">
  <em>⚠️ For research use only. Results were obtained by procedures that were not CLIA validated.</em>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/Pipeline-Sanibel%202.0-blue?style=plastic" />
  <img src="https://img.shields.io/badge/Nextflow-≥22.10-brightgreen?style=plastic&logo=nextflow" />
  <img src="https://img.shields.io/badge/Python-3.10+-yellow?style=plastic&logo=python" />
  <img src="https://img.shields.io/badge/License-Apache%202.0-red?style=plastic" />
</p>

## 🦠🧬 Overview

Sanibel is FL-BPHL's Nextflow bacterial whole-genome sequencing (WGS) analysis pipeline. It performs quality control, *de novo* assembly, species identification, sequence typing and antimicrobial resistance (AMR) detection on paired-end Illumina short reads. Species-specific typing modules run automatically based on species identification results: *Legionella pneumophila* (Legsta), *Klebsiella* (Kleborate), *Shigella* (ShigaTyper), *Streptococcus pyogenes/dysgalactiae* (EMM typing), *Salmonella* (SeqSero2), *E. coli* (SerotypeFinder), *Streptococcus pneumoniae* (SeroBA), *Pseudomonas aeruginosa* (pasty), *Acinetobacter baumannii* (Kaptive), *Vibrio parahaemolyticus* (Kaptive), and *Neisseria meningitidis*/*Haemophilus influenzae* (PMGA). BMGAP2 provides enhanced AMR and antigen analysis for *Neisseria meningitidis* and *Haemophilus influenzae*. PlasmidFinder runs on all samples.


## ⚙️ Dependencies

- **Nextflow** ≥ 22.10 — [installation guide](https://github.com/nextflow-io/nextflow)
- **Apptainer/Singularity** — [installation guide](https://apptainer.org/docs/user/latest/)
- **SLURM** workload manager (This applies only if HiPerGator is used)
- **Conda** (for the SANIBEL environment)


## 💻 Resource Requirements

The pipeline is designed for HPC environments. The most resource-intensive step is **Unicycler** (*de novo* assembly via SPAdes), which runs once per sample.

| Resource | Recommended | Minimum |
|----------|-------------|---------|
| CPUs | 20 (runs 2 assemblies in parallel) | 4 (1 assembly at a time, slower) |
| RAM | 64 GB | 16 GB |
| Disk (input + output) | ~10 GB per sample | — |
| Disk (Kraken2 database) | ~8.5 GB | — |

> **Running locally with fewer CPUs:** Nextflow will still run but your machine will be oversubscribed during assembly. Each Unicycler job requests 10 CPUs by default, on a machine with fewer cores the OS will time-share threads and assembly will complete more slowly but will not fail. You can lower the `cpus` value for `unicycler` in `nextflow.config` to match your hardware.

**Estimated runtime** (12 samples, 20 CPUs, HPC): ~3–4 hours total, dominated by Unicycler (~25 min/sample, 2 running in parallel).


## 🛠️ Setup

### 1. Create the conda environment

```bash
$ conda create -n SANIBEL -c conda-forge -c bioconda python=3.10 pandas=1.5.3 openpyxl=3.1.5 biopython=1.78 mash=2.3 blast=2.17.0
$ conda activate SANIBEL
```

`pandas`, `openpyxl`, `biopython`, `mash`, and `blast` are required by the BMGAP2 modules (used for Nm/Hi samples).

### 2. Configure params.yaml

Edit `params.yaml` and set paths for your environment:

```yaml
# Input / Output absolute path, no trailing slash "/"
input:  "/full/path/to/fastqs"
output: "/full/path/to/output"

# non-HiPerGator users: set this to your BMGAP2 analysis_scripts directory
bmgap2_db:  "/blue/bphl-florida/share/bmgap2"

# non-HiPerGator users: set this to your Kraken2 database directory
kraken_db:  "/blue/bphl-florida/share/kraken2_databases/k2_standard_8GB_20260226"
```

### 3. Configure sanibel.sh

Set `NXF_APPTAINER_CACHEDIR` to your image cache directory and add your email address for job notifications:

```bash
export NXF_APPTAINER_CACHEDIR=/path/to/apptainer/cache
#SBATCH --mail-user=your@email.gov
```

## BMGAP2 Setup

> **HiPerGator users:** BMGAP2 is already installed and configured on the cluster. The `bmgap2_db` path in `params.yaml` is pre-set. Skip this section entirely.

[BMGAP2](https://github.com/CDCgov/BMGAP2) (*Bacterial Meningitis Genome Analysis Pipeline 2*) runs automatically on every sample. The Python scripts it invokes check the MLST scheme internally and skip any sample that is not *N. meningitidis* or *H. influenzae*.

For non-HiPerGator users, BMGAP2 must be installed before running the pipeline. Its scripts run directly on the host (not inside a container) and are invoked by the three `bmgap2_*` Nextflow modules. Follow the [BMGAP2 installation instructions](https://github.com/CDCgov/BMGAP2) to clone the repository and build all required databases, then set `bmgap2_db` in `params.yaml` to the `analysis_scripts` directory.


## How to Run

Place input FASTQ files in the directory specified by `params.input`. Both naming conventions are supported:

| Convention | Example |
|------------|---------|
| Illumina native | `SAMPLE_S1_L001_R1_001.fastq.gz` |
| Simplified | `SAMPLE_1.fastq.gz` |


## 🐊 HiPerGator Usage
```bash
sbatch sanibel.sh
```

## ⚡ Local Usage
```bash
nextflow run sanibel.nf -profile apptainer -params-file params.yaml
```

## Workflow Diagram

```mermaid
flowchart TD
    A[Paired FASTQ Input] --> B[FastQC]
    B --> C[Trimmomatic]
    C --> D[BBTools]
    D --> E[FastQC2]
    B --> F[MultiQC]
    E --> F

    D --> G[Mash]
    D --> H[Unicycler]
    D --> I[Kraken2]

    H --> J[Quast]
    H --> K[AMRFinder]
    H --> L[MLST]
    H --> M[Prokka]

    G --> O[parse_assembly]
    J --> O
    O --> P[readssum]
    D --> P
    P --> Q[parse_reads]
    Q --> M

    M --> R[parse_typing]
    L --> R
    I --> R

    R --> V[Species-specific analyses]

    V -->|Legionella pneumophila| V1[Legsta]
    V -->|Klebsiella| V2[Kleborate]
    V -->|Shigella| V3[ShigaTyper]
    V -->|Streptococcus pyogenes / dysgalactiae| V4[EMM Typing]
    V -->|Salmonella| V5[SeqSero2]
    V -->|E. coli| V6[SerotypeFinder]
    V -->|All samples| V7[PlasmidFinder]
    V -->|Streptococcus pneumoniae| V8[SeroBA]
    V -->|Pseudomonas aeruginosa| V9[pasty]
    V -->|Acinetobacter baumannii| V10[Kaptive AB]
    V -->|Vibrio parahaemolyticus| V11[Kaptive VP]
    V -->|Neisseria meningitidis / H. influenzae| V12[PMGA]
    V12 --> S[BMGAP2 AMR]
    S --> T[BMGAP2 LocusExtractor]
    T --> U[BMGAP2 BMScan]

    V1 --> W[generate_row]
    V2 --> W
    V3 --> W
    V4 --> W
    V5 --> W
    V6 --> W
    V7 --> W
    V8 --> W
    V9 --> W
    V10 --> W
    V11 --> W
    U --> W

    W --> X[sum_report.txt / nm_sum_report.txt / hi_sum_report.txt]

    style S fill:#9cf,stroke:#333
    style T fill:#9cf,stroke:#333
    style U fill:#9cf,stroke:#333
    style V1 fill:#fef,stroke:#333
    style V2 fill:#fef,stroke:#333
    style V3 fill:#fef,stroke:#333
    style V4 fill:#fef,stroke:#333
    style V5 fill:#fef,stroke:#333
    style V6 fill:#fef,stroke:#333
    style V7 fill:#fef,stroke:#333
    style V8 fill:#fef,stroke:#333
    style V9 fill:#fef,stroke:#333
    style V10 fill:#fef,stroke:#333
    style V11 fill:#fef,stroke:#333
    style V12 fill:#fef,stroke:#333
    style X fill:#f96,stroke:#333,stroke-width:3px
```

### Modules

<small>Sanibel is made possible thanks to the following tools:</small>

<small>

| Module | Tool | Version |
|--------|------|---------|
| `fastqc` / `fastqc2` | FastQC | 0.12.1 |
| `trimmomatic` | Trimmomatic | 0.40 |
| `bbtools` | BBTools | 39.77 |
| `multiqc` | MultiQC | 1.33 |
| `mash` | Mash | 2.3 |
| `unicycler` | Unicycler | 0.5.1 |
| `kraken` | Kraken2 | 2.17.1 |
| `quast` | QUAST | 5.3.0 |
| `readssum` | Lyveset | 2.0.1 |
| `prokka` | Prokka | 1.15.6 |
| `amrfinder` | AMRFinderPlus | 4.2.7 |
| `mlst` | MLST | 2.32.2 |
| `pmga` | PMGA | 3.0.2 |
| `bmgap2_amr` / `bmgap2_locusextractor` / `bmgap2_bmscan` | BMGAP2 | — |
| `legsta` | Legsta | 0.5.1 |
| `kleborate` | Kleborate | 3.2.4 |
| `shigatyper` | ShigaTyper | 2.0.5 |
| `emm_typing` | emm-typing-tool | 0.0.1 |
| `seqsero2` | SeqSero2 | 1.3.2 |
| `serotypefinder` | SerotypeFinder | 2.0.2 |
| `plasmidfinder` | PlasmidFinder | 3.0.3 |
| `seroba` | SeroBA | 2.0.5 |
| `pasty` | pasty | 2.2.1 |
| `kaptive_ab` / `kaptive_vp` | Kaptive | 3.2.0 |

</small>

## 📁 Output

All per-sample results are written to `params.output/<sample_id>/`. Depending on which species are in the run, up to three summary files are written to `params.output/`:

| File | Samples | Cols | Key fields |
|------|---------|------|------------|
| `sum_report.txt` | All | 21 | ID · species (Mash/Kraken) · MLST scheme/ST · serotype · QC metrics (reads, coverage, assembly stats, GC, CDS) |
| `nm_sum_report.txt` | *N. meningitidis* only | 26 | ID · PMGA serogroup · BMGAP2 AMR alleles/phenotypes · vaccine antigen coverage (4CMenB) |
| `hi_sum_report.txt` | *H. influenzae* only | 23 | ID · PMGA capsule type · BMGAP2 AMR alleles/phenotypes |

> For Nm and Hi samples, `serotype` in `sum_report.txt` reflects the PMGA prediction. Species-specific details are in the dedicated reports.


## 🤝 Contributing
We welcome contributions to make Sanibel better! Feel free to open issues or submit pull requests to suggest any additional features or enhancements!

## 📧 Contact
**Email**: bphl-sebioinformatics@flhealth.gov

## ⚖️ License
Sanibel is licensed under the [Apache License ](https://github.com/BPHL-Molecular/Sanibel/blob/main/LICENSE).
