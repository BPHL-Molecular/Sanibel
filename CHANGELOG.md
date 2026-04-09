# Changelog

All notable changes to Sanibel are documented in this file.

---

## [2.0.0] — 2026-04-09

### Complete DSL2 Rewrite
The pipeline has been fully rewritten from the ground up in modern Nextflow DSL2.

- Main workflow file renamed from `flaq_amr_plus2.nf` → `sanibel.nf`.
- Input handling moved from manual string enumeration (`Channel.fromList`) to `channel.fromFilePairs`, supporting both Illumina native and simplified file naming in a single workflow.
- All processes now follow the `tuple val(meta), path(...)` pattern with a `meta` map (replaces raw sample name strings), enabling per-sample metadata propagation throughout the pipeline.
- All modules use named outputs (`emit:`) for explicit channel wiring.

### Renamed Modules
The legacy `pyTask*` and `plusAnalyses` modules have been replaced with purpose-named modules:

| 1.3.0 | 2.0.0 |
|-------|-------|
| `pyTask1` | `parse_assembly` |
| `pyTask2` | `parse_reads` |
| `pyTask3` | `parse_typing` |
| `pyTask4` | `generate_row` |
| `plusAnalyses` | `summary_report` + individual species modules (see below) |

### Extracted Species Modules
Species-specific analyses previously bundled as inline bash inside `plusAnalyses.nf` are now independent, containerized Nextflow modules:

- `legsta` — *Legionella pneumophila* typing (Legsta 0.5.1)
- `kleborate` — *Klebsiella* screening (Kleborate 3.2.4)
- `shigatyper` — *Shigella* serotyping (ShigaTyper 2.0.5)
- `emm_typing` — Group A *Streptococcus* emm typing (emm-typing-tool 0.0.1)
- `seqsero2` — *Salmonella* serotyping (SeqSero2 1.3.2)
- `serotypefinder` — *E. coli* serotyping (SerotypeFinder 2.0.2)
- `plasmidfinder` — Plasmid detection (PlasmidFinder 3.0.3)

### New Modules
Four new species-specific typing modules added:

- **`seroba`** — *Streptococcus pneumoniae* serotyping (SeroBA 2.0.5)
- **`pasty`** — *Pseudomonas aeruginosa* serogroup typing (pasty 2.2.1)
- **`kaptive_ab`** — *Acinetobacter baumannii* K/OC locus typing (Kaptive 3.2.0)
- **`kaptive_vp`** — *Vibrio parahaemolyticus* K/O locus typing (Kaptive 3.2.0)

### Re-enabled Modules
- **`pmga`** — Was present in 1.3.0 but commented out of the workflow. Now runs as a dedicated module (PMGA 3.0.2) for *Neisseria meningitidis* and *H. influenzae* samples.

### BMGAP2 Integration
Three new modules added for enhanced *N. meningitidis* and *H. influenzae* analysis:

- **`bmgap2_amr`** — Mutation-based AMR profiling
- **`bmgap2_locusextractor`** — Vaccine antigen identification
- **`bmgap2_bmscan`** — Species confirmation

BMGAP2 runs automatically on every sample; the Python scripts check the MLST scheme internally and skip non-Nm/Hi samples. No `meningitis` flag is required. `params.bmgap2_db` defaults to `/blue/bphl-florida/share/bmgap2` (HiPerGator path).

### Container Version Updates
| Module | 1.3.0 | 2.0.0 |
|--------|-------|-------|
| `fastqc` / `fastqc2` | 0.11.9 | 0.12.1 |
| `trimmomatic` | 0.39 | 0.40 |
| `bbtools` | 38.76 | 39.77 |
| `multiqc` | 1.8 | 1.33 |
| `mash` | 2.2 | 2.3 |
| `unicycler` | 0.4.7 | 0.5.1 |
| `kraken` | 2.0.8-beta | 2.17.1 |
| `quast` | 5.0.2 | 5.3.0 |
| `readssum` (lyveset) | 1.1.4f | 2.0.1 |
| `prokka` | 1.14.5 | 1.15.6 |
| `amrfinder` | 3.10.1 | 4.2.7 |
| `mlst` | 2.19.0 | 2.32.2 |
| `pmga` | latest *(unused)* | 3.0.2 |
| `kleborate` | 2.2.0 | 3.2.4 |
| `shigatyper` | 2.0.1 | 2.0.5 |
| `seqsero2` | 1.2.1 | 1.3.2 |
| `serotypefinder` | 2.0.1 | 2.0.2 |
| `plasmidfinder` | 2.1.6 | 3.0.3 |

### Breaking Changes
- **`kleborate`** updated for v3 CLI: new flags (`-a`, `-o kleborate_out`, `-p kpsc`, `--trim_headers`); output path changed to `kleborate_out/klebsiella_pneumo_complex_output.txt`.
- **`plasmidfinder`** updated for v3: entry point changed from `plasmidfinder.py` to `python -m plasmidfinder`.

### Configuration Changes
- Single `nextflow.config` replaces the `configs/config_template.config` split used in 1.3.0.
- Added `profiles` block supporting `standard`, `docker`, `singularity`, and `apptainer` execution profiles.
- `autoMounts = true` moved inside the `singularity {}` and `apptainer {}` profile blocks (was incorrectly at top level / in docker block in 1.3.0).
- CPU and memory resource limits added for every process.
- `params.bmgap2_db` default added.

### Bug Fixes
- Fixed `prefix` variable scope error in 10 modules (`amrfinder`, `trimmomatic`, `bbtools`, `mash`, `unicycler`, `kraken`, `mlst`, `pmga`, `readssum`, `prokka`): `${prefix}` in `output:` blocks replaced with `${meta.id}`, which is in scope.
- Fixed 19 Nextflow DSL2 linter errors in `sanibel.nf`: moved top-level statements inside `workflow {}`, renamed `Channel` → `channel` (lowercase), and prefixed unused closure parameters with `_`.

---

## [1.3.0] — 2025-10-16

Initial tagged release.
