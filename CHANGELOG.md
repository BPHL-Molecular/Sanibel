# Changelog

All notable changes to Sanibel are documented in this file.

---

## [Unreleased]

### Species ID QC — Confidence Gating + QC Columns
Added quality-control gating and reporting for species identification and assembly quality.

- **`modules/skani.nf`** — Species routing now requires the top hit to clear both ANI (`skani_routing_min_ani`, default 95) **and** query alignment fraction (`skani_routing_min_af`, default 50). The AF gate rejects high-identity hits that cover only part of the genome.
- **`sanibel.nf`** — When skani returns no confident species, `meta.species` / `meta.genus` are set to `Unknown` instead of falling back to the Mash top hit, so species-specific typing (serotypers, PMGA, BMGAP2) is withheld for those samples. MLST, AMRFinder, and PlasmidFinder are unaffected. Prevents mis-routing on divergent/mixed samples (e.g. a Mash top hit of *Staphylococcus* on a non-ID *Bacillaceae*).
- **`bin/summary_report.py`** — Two new `sum_report.txt` columns:
  - `species_id_qc` (after `skani_reference`): `Pass` / `No WGS ID (ANI <95%)` / `No WGS ID (AF <50%)` / `NO ID (ANI <80%)`, derived from the top skani hit's ANI and alignment fraction. `skani_species` still shows the raw top hit for review.
  - `assembly_qc` (last column): `Pass` / `Warn` / `Fail` with reasons, from coverage (`<40x` Fail), contig count (`>=200` Warn, `>500` Fail), and N50 (`<15000` Warn).
- **`nextflow.config`** / **`modules/summary_report.nf`** — New tunable params `skani_routing_min_af` (50), `qc_min_coverage` (40), `qc_warn_contigs` (200), `qc_fail_contigs` (500), `qc_min_n50` (15000), passed through to `summary_report.py`.

### Maintenance — Behavior-Preserving Refactor
Internal cleanup to reduce duplication. No change to pipeline outputs.

- **`bin/sanibel_taxonomy.py`** — Now the single home for the shared 16S-BLAST and Kraken row parsers (`iter_blast16s_rows`, `iter_kraken_species_rows`), the 16S thresholds (`BLAST16S_MIN_LENGTH`, `BLAST16S_MIN_PIDENT`), the accession regex (`ACCESSION_RE` / `find_accession`), and `species_of`. `aggregate_species_id.py`, `build_candidate_pool.py`, `parse_assembly.py`, and `summary_report.py` import these instead of each carrying their own copy.
- **`bin/build_candidate_pool.py`** — Dropped the unused `sample_id` argument; **`modules/build_candidates.nf`** no longer passes it.
- **`bin/summary_report.py`** — Collapsed the three duplicated standard-row blocks (neisseria / hinfluenzae / other) into one.
- **`bin/bmgap2_helpers.py`** (new) — Shared MLST-scheme / sample-name boilerplate for the three `run_bmgap2_*.py` host scripts. The per-script meningitis re-check is preserved.
- **`modules/kaptive.nf`** (new) — Single parameterized Kaptive module (`variant` of `ab` / `vp`), invoked via include aliases; replaces `kaptive_ab.nf` and `kaptive_vp.nf`. Output filenames and publish dirs unchanged.
- **`sanibel.nf`** — Added a `rebind()` helper for the repeated re-key / join / re-emit idiom.
- **`nextflow.config`** — Added process-level `cpus` / `memory` defaults (per-process blocks override where they differ) and a `withName: 'bmgap2_.*'` selector for the shared BMGAP2 settings. Removed the redundant module-level `errorStrategy` from `candidate_references.nf`.

### Output — Smaller per-sample output
- **`modules/candidate_references.nf`** — No longer publishes the downloaded RefSeq reference genomes to `<sample>/candidate_references/`. They stay internal to the run (skani still reads them from `work/`), removing the largest per-run output artifact. Reference provenance remains in `_skani.tsv`.
- **`modules/skani.nf`** — Publishes only `_skani.tsv`; the internal `_skani_species.txt` routing file is no longer copied to `<sample>/skani/`. Its content still drives species routing.

### Species Identification — Candidate Pool + skani ANI Confirmation
Replaced the previous Kraken2/Mash agreement heuristic with a two-stage workflow:
three tools (Mash and Kraken2 on the reads, 16S rRNA BLAST on the assembly) nominate
a ranked pool of candidate species, multiple RefSeq reference genomes (N=5 by
default) are downloaded per candidate, and `skani` confirms the species by ANI. New `sum_report.txt` columns:
`blast_16s_tophit`, `blast_16s_pident`, `skani_species`, `skani_ani`,
`skani_align_fraction`, `skani_reference`, `contamination_flag`.

- **`modules/blast_16s.nf`** — Downloads the NCBI `16S_ribosomal_RNA` BLAST database once (cached via `storeDir`) and megablasts each assembly. `perc_identity` set to 97 to match downstream filters.
- **`bin/build_candidate_pool.py`** / **`modules/build_candidates.nf`** — Builds the ranked candidate pool from all three tools. Includes every distinct Mash species and every qualifying 16S species; Kraken2 candidates now require a minimum clade-read count (`KRAKEN_MIN_READS = 10`) so low-abundance `0.00%` taxa no longer enter the pool. Pool capped at 15.
- **`modules/candidate_references.nf`** — Downloads one RefSeq genome per candidate via `datasets` (accession first, taxon-name fallback). Downloads parallelized (`xargs -P 3`) with `maxForks 4` to bound concurrent NCBI sessions; `errorStrategy 'ignore'`.
- **`modules/skani.nf`** — ANI confirmation against all downloaded references at once; the best ANI hit drives the species call.
- **`bin/aggregate_species_id.py`** / **`modules/aggregate_species_id.nf`** — Retained to supply `contamination_flag` (a foreign 16S genus on overlapping contigs), with a synonymous-genus table to avoid false positives between 16S-indistinguishable genera. Output published to `candidate_species/`.
- **`bin/summary_report.py`** — Reports 16S at genus level as `Genus spp.` (`blast_16s_tophit`), cross-checked against Mash and Kraken so a contaminant contig is not reported as the top hit. Removed the unused `skani_notes` helper.
- **`modules/mash.nf`** — Widened the distance output to the top 50 hits (deduplicated by species downstream) and renamed the artifact to `*_mash_distances.tab`.

### Mash Reference Sketch Refresh
- **`nextflow.config`** — Mash container pinned to `staphb/mash:2.3-RefSeqProkv235`, baking in a current (2026) RefSeq prokaryote sketch from `update_mash_dist`, replacing the stale 2019 gembox sketch.
- **`bin/parse_assembly.py`** / **`bin/build_candidate_pool.py`** — Parse the new `Genus_species_<ACCESSION>` sketch naming and recover the accession via the existing `_ACC_RE` regex (legacy `-.-` format still supported).

### skani as Species-ID and Contamination Arbiter
Made skani the source of truth for the reported organism, the 16S anchor, and the contamination flag. The report previously inherited these from Mash, which is unstable on mixed or divergent samples (e.g. a contaminated *Acinetobacter* flipping to *Staphylococcus* and inverting the contamination flag).

- **`bin/sanibel_taxonomy.py`** (new): single home for the 16S-synonymous genus table and the contig-overlap contamination logic, imported by both `aggregate_species_id.py` and `summary_report.py` so the two cannot drift.
- **`bin/aggregate_species_id.py`**: imports the shared module; contamination detection now also suppresses synonymous-genus pairs (e.g. *Escherichia*/*Shigella*), matching the table's stated intent.
- **`bin/summary_report.py`**: `blast_16s_tophit` is anchored on the skani genus; `contamination_flag` is recomputed against the skani genus; the reported organism is relabeled to the Mash/Kraken consensus when skani's winner is a 16S-synonymous genus of it (*Escherichia*/*Shigella*). Each path falls back to the prior Mash/Kraken behavior when skani returns no result.

### Multi-Reference skani ANI
Recovers ANI above the 95% species boundary for clinical isolates whose representative RefSeq strain is divergent (*E. coli* vs K-12, *L. monocytogenes* vs EGD-e). Supersedes the single-genome download described above.

- **`modules/candidate_references.nf`**: downloads up to N RefSeq genomes per candidate species instead of one. Lists accessions with `datasets summary genome taxon ... --limit N` (`--limit` is not valid on the taxon download), always includes the sketch representative, fetches them in a single `datasets download genome accession` call, and names each file `Genus_species__<accession>.fna`. skani picks the best ANI across all strains; keeping the representative makes the reference set a strict superset, so ANI cannot regress versus the single-genome behavior.
- **`bin/summary_report.py`**: `parse_skani` reads the new per-accession filenames, reporting a clean `Genus_species` label and the winning strain's accession.
- **`nextflow.config`**: new `params.refseq_refs_per_candidate` (default 5).

### Added
- **`modules/lissero.nf`** — New LisSero serogroup-typing module for *Listeria monocytogenes*. Runs on the assembly and is gated by `meta.species == 'Listeria_monocytogenes'`, mirroring the other species-specific modules. Container: `staphb/lissero:0.4.10`.
- **`sanibel.nf`** — Included and invoked `lissero`; added `lissero.out.done` to the summary-report barrier.
- **`nextflow.config`** — Added `withName: lissero` resource/container block (`errorStrategy = 'ignore'`).
- **`summary_report.py`** — Added `get_listeria_serotype()` parser; the LisSero `SEROTYPE` value now populates the `serotype` column of `sum_report.txt`.

### Configuration
- **`nextflow.config`** — `bmgap2_amr`, `bmgap2_locusextractor`, and `bmgap2_bmscan` set `cache = false` to force re-execution on every run.

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
