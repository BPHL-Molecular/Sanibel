# Changelog

All notable changes to Sanibel are documented in this file.

---

## [2.0.0] — 2026-04-22

### Bug Fixes
- **`kaptive_ab`** — Fixed command for Kaptive v3: executable renamed from `kaptive.py` to `kaptive`; replaced file-path DB arguments with built-in keywords `ab_k` and `ab_o`.
- **`emm_typing`** — Fixed EMBOSS library conflict preventing execution. Added `containerOptions "--bind ${task.workDir}:/EMBOSS-6.6.0/emboss/.libs"`.
- **`sanibel.nf`** — Fixed `_` reserved identifier in two closures (renamed to `_id` / `_ids`).
- **`sanibel.nf`** — Fixed `mash_species` Groovy GString assignment that caused all species-specific `.filter {}` comparisons to return `false`. Changed to plain String concatenation (`fields[0] + '_' + fields[1]`). Affected modules: `emm_typing`, `seroba`, `pasty`, `kaptive_ab`, `kaptive_vp`, `serotypefinder`.
- **`summary_report.py`** — MLST `-` values now normalised to `Not detected`.
- **`summary_report.py`** — LocusExtractor CSV parser now prefers prokka-annotated rows for complete antigen ORF data.
- **`summary_report.py`** — `normalize_le_value()` extended to handle `Incomplete ORF`, `New-BLASTonly`, and `New-PCR` result strings.
- **`summary_report.py`** — Added missing serotype parsers: `get_pneumococcal_serotype()` (SeroBA), `get_acinetobacter_serotype()` (Kaptive), `get_pseudomonas_serotype()` (pasty).
- **`summary_report.py`** — Report column order updated: species ID, MLST, and serotype now precede QC metrics in `sum_report.txt`. Nm/Hi species-specific reports (`nm_sum_report.txt`, `hi_sum_report.txt`) condensed to species-only columns.
- **`params.yaml`** — Added `kraken_db` parameter.

### Complete DSL2 Rewrite
The pipeline has been fully rewritten in modern Nextflow DSL2.

- Main workflow file renamed from `flaq_amr_plus2.nf` → `sanibel.nf`.
- Input handling moved from manual string enumeration (`Channel.fromList`) to `channel.fromFilePairs`, supporting both Illumina native and simplified file naming in a single workflow.
- All processes now follow the `tuple val(meta), path(...)` pattern with a `meta` map (replaces raw sample name strings), enabling per-sample metadata propagation throughout the pipeline.
- All modules use named outputs (`emit:`) for explicit channel wiring.

### Renamed / Replaced Modules
The legacy `pyTask*` and `plusAnalyses` modules have been replaced:

| 1.3.0 | 2.0.0 |
|-------|-------|
| `pyTask1` | `parse_assembly` |
| `pyTask2` | *(removed — read metrics now sourced via `readssum` + `meta.genome_size`)* |
| `pyTask3` | *(removed — typing data threaded via individual files)* |
| `pyTask4` | *(removed — replaced by `summary_report`)* |
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

### BMGAP2 Integration
Three new modules added for enhanced *N. meningitidis* and *H. influenzae* analysis:

- **`bmgap2_amr`** — Mutation-based AMR profiling
- **`bmgap2_locusextractor`** — Vaccine antigen identification
- **`bmgap2_bmscan`** — Species confirmation

BMGAP2 runs automatically on every sample; the python scripts check the MLST scheme internally and skip non-Nm/Hi samples. `params.bmgap2_db` defaults to `/blue/bphl-florida/share/bmgap2` (HiPerGator path).

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
| `pmga` | latest | 3.0.2 |
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
- CPU and memory resource added for every process.
- `params.bmgap2_db` default added.

### Bug Fixes
- Fixed 19 Nextflow DSL2 linter errors in `sanibel.nf`: moved top-level statements inside `workflow {}`, renamed `Channel` → `channel` (lowercase), and prefixed unused closure parameters with `_`.
- `bbtools_phix`: fixed `-in`/`-in2` flag order causing read-pair swap.
- `fastqc2`: fixed input declared as `path` instead of `tuple val(meta), path`.
- `kraken`: fixed container path and database mount.
- `prokka`: added `export _JAVA_OPTIONS="-XX:-UsePerfData"` to suppress JVM crash on HiPerGator.
- `unicycler`: reduced memory request to avoid SLURM OOM kills.
- `readssum`: now uses `meta.genome_size` directly instead of reading pyoutputs.
- `sanibel.nf`: fixed closure parameter warnings — unused params prefixed with `_`; params shadowing imported process names renamed.

### Architecture — Eliminated pyoutputs Accumulation Pattern
The in-place CSV accumulation pattern (`parse_assembly` → `parse_reads` → `parse_typing` mutating a shared file per sample) has been fully removed and replaced with a clean Nextflow data-flow approach.

#### Meta Enrichment
After `parse_assembly`, the workflow now enriches `meta` inline:
- `meta.mash_genus` — Mash top-hit genus (e.g. `Salmonella`)
- `meta.mash_species` — Mash top-hit genus + species (e.g. `Salmonella_enterica`)
- `meta.genome_size` — total assembly length (used by `readssum` for coverage calculation)

#### Summary Report Rewrite
`bin/summary_report.py` completely rewritten:
- Discovers sample IDs from `*_assembly_stats.txt` files staged in the work directory
- Parses each mandatory staged file directly (`_assembly_stats.txt`, `_readMetrics.txt`, prokka `.txt`, `.mlst`, Kraken `.report`, `sta.txt`)
- Reads species-specific outputs from `--outdir/{sample_id}/tool/` (serotypefinder, kleborate, seqsero2, emm_typing, shigatyper, legsta)
- BMGAP2: checks for `bmgap2_amr/{sample_id}*amr_data.json` — if present, parses all three BMGAP2 output dirs
- Routes by MLST scheme: `sum_report.txt` (standard), `sum_report_nm.txt` (neisseria), `sum_report_hi.txt` (hinfluenzae)

#### Species Gating
All 11 species-specific modules now gate via `meta.mash_genus` / `meta.mash_species` Groovy interpolation in their bash blocks. Inputs no longer include a pyoutputs path.

#### BMGAP2 Updates
- All three BMGAP2 modules now thread `mlst_file` through the chain (was pyoutputs)
- Python scripts now read MLST scheme from the `.mlst` file (whitespace field[1]) instead of a pyoutputs CSV field
- Fixed `assembly_dir` path: `{id}_assembly/` → `assembly/`
- `params.meningitis` flag removed — BMGAP2 self-gates via MLST scheme

### Output Directory Standardization
| Output | Before | After |
|--------|--------|-------|
| Unicycler assembly | `{sample_id}_assembly/` | `assembly/` |
| AMRFinder | *(inline)* | `amrfinder/` |
| MLST | *(inline)* | `mlst/` |
| SeqSero2 | `seqsero2_output/` | `seqsero2/` |
| Kleborate | `kleborate_output/` | `kleborate/` |
| SerotypeFinder | `serotypefinder_output/` | `serotypefinder/` |
| emm-typing | `emm_output/` | `emm_typing/` |

---

## [1.3.0] — 2025-10-16

Initial tagged release.
