# Changelog

All notable changes to Sanibel are documented in this file.

---

## [Unreleased]

### Nextflow 26 support
- `modules/*.nf`: dynamic `publishDir` directives rewritten as closures for the v2 strict parser.
- `sanibel.sh`: loads the default `nextflow` module instead of pinning `nextflow/25.10.4`.
- Parameters are set via `-params-file`.

### Species identification
Candidate pool plus skani ANI replaces the Kraken2/Mash agreement heuristic.

- `modules/blast_16s.nf`: NCBI `16S_ribosomal_RNA` megablast per assembly, `perc_identity` 97.
- `bin/build_candidate_pool.py`: ranked pool, each tool's top 3 seeded, capped at 15; Kraken2 needs 10+ clade reads.
- `modules/candidate_references.nf`: 5 RefSeq genomes per candidate plus the Mash sketch representative.
- `modules/skani.nf`: species call requires ANI >= 95 and alignment fraction >= 50.
- `sanibel.nf`: no confident skani call sets `meta.species` to `Unknown` and withholds species-specific typing.
- `bin/summary_report.py`: `skani_species` is skani's best reference; `blast_16s_tophit` anchored on the skani genus.
- `bin/sanibel_taxonomy.py` (new): shared synonymous-genus table and contig-overlap contamination logic.
- Contamination reported once in `sum_report.txt`, against the skani genus; `SYNONYMOUS_PAIRS` suppresses 16S-indistinguishable pairs.
- `candidate_species.txt`: 2-of-3 vote review record; does not feed the summary report.
- `nextflow.config`: Mash pinned to `staphb/mash:2.3-RefSeqProkv235` (2026 RefSeq sketch).
- `modules/mash.nf`: top 50 distances, artifact renamed `*_mash_distances.tab`.

### Report columns
- New: `blast_16s_tophit`, `blast_16s_pident`, `skani_species`, `skani_ani`, `skani_align_fraction`, `skani_reference`, `contamination_flag`.
- `species_id_qc` (col 2): `PASS` / `REVIEW (borderline ANI)` / `NO ID (ANI < 95%)`.
- `assembly_qc` (col 4): `PASS` / `Warning:` / `FAIL:` / `FAIL (Contamination)`, from coverage, contigs, N50 and contamination.
- Kraken2 columns renamed `kraken2_species` / `kraken2_percent`; output moved to `<sample>/kraken2/`.
- AMR columns renamed `amr_gene_symbol` / `amr_subclass`.
- QC thresholds are hardcoded constants, not `nextflow.config` params.

### Reliability
- `sanibel.nf`: the optional-typing barrier defaults to `true`, so a run with no species-specific output still produces a report.

### BMGAP2 resume
- `modules/bmgap2_amr.nf`: takes the PMGA BLAST JSON as a staged `path` input instead of reading `${params.output}/<id>/pmga`.
- `modules/bmgap2_locusextractor.nf`, `modules/bmgap2_bmscan.nf`: take the contigs and the prokka annotation as staged `path` inputs and build the scan directory in the work directory. Previously they read `${params.output}/<id>/assembly` with no DAG edge to prokka, so LocusExtractor could run before its preferred prokka FASTA existed.
- All three declare their results as process outputs and publish them; `cache = false` dropped from the `withName: 'bmgap2_.*'` selector, so a resume republishes them from `work/`.
- `modules/pmga.nf`: second output named `emit: files`.
- `bin/run_bmgap2_*.py`: explicit `--sample` / `--mlst` / `--outdir` / `--db` arguments; no reference to `params.output`. `run_bmgap2_amr.py` no longer writes a `json/` subdirectory into PMGA's published output.
- `bin/bmgap2_helpers.py`: `read_scheme_and_sample()` replaced by `read_scheme()`, which takes the sample ID rather than deriving it from the output path.
- `bin/run_bmgap2_bmscan.py`: BMScan's JSON renamed to `<sample>_species_analysis.json` so many samples can be staged into one `summary_report` task.
- `modules/summary_report.nf`, `bin/summary_report.py`: BMGAP2 results arrive as staged inputs; `parse_bmgap2()` reads the work directory instead of `${params.output}/<id>`.

### Output
- `modules/candidate_references.nf`: reference genomes no longer published; accession manifest published instead.
- `modules/skani.nf`: publishes only `_skani.tsv`.

### Added
- `modules/lissero.nf`: LisSero serogroup typing for *Listeria monocytogenes*, populating the `serotype` column.

### Maintenance
- `bin/sanibel_taxonomy.py`: single home for the shared 16S/Kraken parsers, thresholds, accession regex and `genus_of`.
- `bin/bmgap2_helpers.py` (new): shared boilerplate for the three `run_bmgap2_*.py` scripts.
- `modules/kaptive.nf` (new): parameterized module replacing `kaptive_ab.nf` / `kaptive_vp.nf`.
- `bin/summary_report.py`: collapsed the three duplicated standard-row blocks into one.
- `bin/build_candidate_pool.py`: dropped the unused `sample_id` argument.
- `sanibel.nf`: added a `rebind()` helper for the repeated re-key / join idiom.
- `nextflow.config`: process-level `cpus` / `memory` defaults; `withName: 'bmgap2_.*'` selector; `bmgap2_*` set `cache = false`.

---

## [2.0.1] — 2026-06-19

### Improvements
- **`environment.yaml`** — Added conda environment file for reproducible setup via `conda env create -f environment.yaml`; replaces the manual package-list install command.
- **`sanibel.sh`** — Pinned Nextflow module to `nextflow/25.10.4`; Nextflow ≥ 26.0 is not supported due to DSL2 breaking changes.
- **`sanibel.sh`** — Removed `./cache` from the cleanup line; `NXF_APPTAINER_CACHEDIR` is now the single persistent cache location.

### README
- Added Nextflow version constraint warning: ≥ 26.0 not supported.
- Added git clone step to setup instructions.
- Added Kraken2 and BMGAP2 setup section for non-Florida-BPHL users.
- Generalized "HiPerGator" language to "Florida BPHL" throughout.
- Changed mermaid flowchart direction from left-right to top-down.

---

## [2.0.0] — 2026-05-04

### Bug Fixes
- **`bin/parse_assembly.py`** — Fixed parsing of MASH reference IDs containing two `-.-` separators (e.g. some *Listeria monocytogenes* RefSeq entries: `...-PRJNA224116-.-GCF_000307085.1-.-Listeria_monocytogenes_...`). The previous `split('-.-', 1)` produced a wrong accession (`PRJNA224116`) and a malformed genus (`GCF`). The fix splits on all `-.-` occurrences and selects `segs[-2]`/`segs[-1]` as accession/organism for entries with ≥ 3 segments.
- **`bin/parse_assembly.py`** — Fixed genus/species extraction for reference organisms whose names begin with an underscore (e.g. `_Cellvibrio_japonicus`). Empty strings produced by leading underscores are now filtered out before indexing the parts list.
- **`summary_report.py`** — Added missing `get_vibrio_serotype()` implementation; *Vibrio parahaemolyticus* serotype from Kaptive VP was previously never reported in `sum_report.txt`.
- **`summary_report.py`** — Fixed NM/HI summary report column order: `bmgap2_species`, `bmgap2_mlst_st`, `bmgap2_mlst_cc`, and `serotype_notes` now precede `nm_serogroup`/`hi_serotype` and `predicted_resistance`.
- **`bin/run_bmgap2_locusextractor.py`** — Fixed `PermissionError` on `allele_reference.txt` caused by file ownership after BMGAP2 refactor; resolved on HiPerGator with `chmod g+rw`.
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

### Improvements
- **`modules/pmga.nf`** — Added `--threads ${task.cpus}` flag; PMGA previously defaulted to 1 thread regardless of the CPU allocation in `nextflow.config`.
- **`nextflow.config`** — Increased BMGAP2 process CPU/memory allocations (`bmgap2_amr`: 4 CPU / 4 GB; `bmgap2_locusextractor`: 4 CPU / 8 GB; `bmgap2_bmscan`: 4 CPU / 4 GB). Added explicit resource block for `pmga` (4 CPU / 8 GB).
- **`sanibel.sh`** — Increased SLURM allocation to 40 CPUs / 200 GB, enabling 4 Unicycler assemblies to run in parallel (down from ~3–4 hrs to ~75 min for a 12-sample run). Added post-run timestamp rename of the output directory.
- **`nextflow.config`** — Default `kraken_db` reverted to `minikraken2_v1_8GB_201904`; the newer `k2_standard_8GB_20260226` database showed lower species-level classification rates due to increased k-mer ambiguity from a larger reference genome set compressed into the same 8 GB space.

### README
- Reorganised **Modules** section: Lyveset moved to Quality Control; PMGA and BMGAP2 moved to Species-Specific; added new **AMR & Mobile Genetic Elements** category for AMRFinderPlus and PlasmidFinder.
- Clarified disk requirement note to specify it applies to the minikraken2 / standard-8 databases.

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
| `bbtools` | 38.76 | 39.34 |
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
