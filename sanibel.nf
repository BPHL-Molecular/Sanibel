#!/usr/bin/env nextflow

/*
  Sanibel Pipeline (named after Sanibel Island in southwest Florida)
  Florida's BPHL Nextflow pipeline for Bacterial WGS Analysis
  Authors: Sarah Schemedes, Yibo Dong, Arnold Rodriguez-Hilario, Molly Mitchell
  Email: bphl-sebioinformatics@flhealth.gov 
*/

nextflow.enable.dsl = 2

include { fastqc }                from './modules/fastqc.nf'
include { trimmomatic }           from './modules/trimmomatic.nf'
include { bbtools_adapters }      from './modules/bbtools.nf'
include { bbtools_phix }          from './modules/bbtools.nf'
include { fastqc2 }               from './modules/fastqc2.nf'
include { multiqc }               from './modules/multiqc.nf'
include { mash }                  from './modules/mash.nf'
include { unicycler }             from './modules/unicycler.nf'
include { kraken }                from './modules/kraken.nf'
include { quast }                 from './modules/quast.nf'
include { parse_assembly }        from './modules/parse_assembly.nf'
include { readssum }              from './modules/readssum.nf'
include { prokka }                from './modules/prokka.nf'
include { amrfinder }             from './modules/amrfinder.nf'
include { mlst }                  from './modules/mlst.nf'
include { pmga }                  from './modules/pmga.nf'
include { download_16s_db }       from './modules/blast_16s.nf'
include { blast_16s }             from './modules/blast_16s.nf'
include { aggregate_species_id }  from './modules/aggregate_species_id.nf'
include { build_candidates }      from './modules/build_candidates.nf'
include { refseq_references }     from './modules/refseq_references.nf'
include { skani }                 from './modules/skani.nf'
include { bmgap2_amr }            from './modules/bmgap2_amr.nf'
include { bmgap2_locusextractor } from './modules/bmgap2_locusextractor.nf'
include { bmgap2_bmscan }         from './modules/bmgap2_bmscan.nf'
include { legsta }                from './modules/legsta.nf'
include { kleborate }             from './modules/kleborate.nf'
include { shigatyper }            from './modules/shigatyper.nf'
include { emm_typing }            from './modules/emm_typing.nf'
include { seqsero2 }              from './modules/seqsero2.nf'
include { serotypefinder }        from './modules/serotypefinder.nf'
include { plasmidfinder }         from './modules/plasmidfinder.nf'
include { seroba }                from './modules/seroba.nf'
include { pasty }                 from './modules/pasty.nf'
include { kaptive_ab }            from './modules/kaptive_ab.nf'
include { kaptive_vp }            from './modules/kaptive_vp.nf'
include { lissero }               from './modules/lissero.nf'
include { summary_report }        from './modules/summary_report.nf'

workflow {
    log.info """
    Sanibel — Bacterial WGS Analysis Pipeline
    ==========================================================================
    input dir   : ${params.input}
    output dir  : ${params.output}
    bmgap2 db   : ${params.bmgap2_db}
    kraken db   : ${params.kraken_db}
    ==========================================================================
    """

    // FASTQ Input channel
    ch_reads = channel.fromFilePairs(
        ["${params.input}/*_{1,2}.fastq.gz",
         "${params.input}/*_R{1,2}_*.fastq.gz"],
        checkIfExists: false
    )
    .map { id, files ->
        def clean_id = id.replaceAll(/_S\d+_L\d+$/, '')
        def meta = [ id: clean_id, single_end: false ]
        [ meta, files ]
    }

    // MLST CC reference tables
    def mlstTablesDir = "${projectDir}/db/mlst_tables"
    [
        [ file("${mlstTablesDir}/neisseria.txt"),
          "https://raw.githubusercontent.com/tseemann/mlst/master/db/pubmlst/neisseria/neisseria.txt" ],
        [ file("${mlstTablesDir}/hinfluenzae.txt"),
          "https://raw.githubusercontent.com/tseemann/mlst/master/db/pubmlst/hinfluenzae/hinfluenzae.txt" ]
    ].each { dst, src ->
        if ( !dst.exists() ) {
            log.info "Downloading ${dst.name} from PubMLST..."
            dst.text = src.toURL().text
        }
    }

    ch_neisseria_txt   = channel.value(file("${mlstTablesDir}/neisseria.txt",   checkIfExists: true))
    ch_hinfluenzae_txt = channel.value(file("${mlstTablesDir}/hinfluenzae.txt", checkIfExists: true))

    // QC & read preprocessing
    ch_fastqc  = fastqc(ch_reads)
    ch_trimmed = trimmomatic(ch_reads)
    ch_adaptertrim = bbtools_adapters(ch_trimmed.reads)
    ch_clean       = bbtools_phix(ch_adaptertrim.reads)
    ch_fastqc2 = fastqc2(ch_clean.reads)

    // Per-sample MultiQC (raw and clean FastQC reports)
    multiqc(
        ch_fastqc.report
            .join(ch_fastqc2.report, by: 0)
            .map { meta, raw_reports, clean_reports ->
                [meta, (raw_reports instanceof List ? raw_reports : [raw_reports]) +
                       (clean_reports instanceof List ? clean_reports : [clean_reports])]
            }
    )

    // Species ID, assembly and read classification
    ch_mash     = mash(ch_clean.reads)
    ch_assembly = unicycler(ch_clean.reads)
    ch_kraken   = kraken(ch_clean.reads)
    ch_quast    = quast(ch_assembly.assembly)

    // Assembly statistics
    ch_parse_assembly = parse_assembly(
        ch_mash.distances
            .join(ch_quast.report, by: 0)
    )

    // Meta + assembly stats channel
    ch_stats = ch_parse_assembly.map { meta, stats ->
        def fields        = stats.text.trim().split(',')
        def enriched_meta = meta + [
            mash_genus:   fields[0],
            mash_species: fields[0] + '_' + fields[1],
            genome_size:  fields[9].toLong()
        ]
        [ enriched_meta, stats ]
    }

    // Meta channel keyed by sample ID
    ch_meta_by_id = ch_stats.map { meta, _stats -> [ meta.id, meta ] }

    // Rebind clean reads and assembly with enriched meta
    ch_clean_enriched = ch_clean.reads
        .map  { meta, reads -> [ meta.id, reads ] }
        .join(ch_meta_by_id)
        .map  { _id, reads, emeta -> [ emeta, reads ] }

    ch_assembly_enriched = ch_assembly.assembly
        .map  { meta, asm -> [ meta.id, asm ] }
        .join(ch_meta_by_id)
        .map  { _id, asm, emeta -> [ emeta, asm ] }

    // Read metrics
    ch_readssum = readssum(ch_clean_enriched)

    // Annotation and typing
    ch_prokka = prokka(
        ch_assembly_enriched
            .map  { meta, asm -> [ meta.id, meta, asm ] }
            .join(ch_stats.map { meta, stats -> [ meta.id, stats ] })
            .map  { _id, meta, asm, stats -> [ meta, asm, stats ] }
    )
    amrfinder(ch_assembly_enriched)
    ch_mlst = mlst(ch_assembly_enriched)
    ch_pmga = pmga(
        ch_assembly_enriched
            .filter { meta, _f -> meta.mash_genus in ['Neisseria', 'Haemophilus'] }
            .join(ch_mlst.out, by: 0)
    )

    // Kraken output with enriched meta
    ch_kraken_enriched = ch_kraken.out
        .map  { meta, report -> [ meta.id, report ] }
        .join(ch_meta_by_id)
        .map  { _id, report, emeta -> [ emeta, report ] }

    // 16S BLAST (DB downloaded once, cached via storeDir)
    ch_16s_db    = download_16s_db().db
    ch_blast_16s = blast_16s(ch_assembly_enriched, ch_16s_db)

    // 2-of-3 species vote (Mash + Kraken2 + 16S BLAST)
    ch_aggregate = aggregate_species_id(
        ch_stats.map { meta, stats -> [ meta.id, meta, stats ] }
            .join(ch_kraken.out.map    { meta, r -> [ meta.id, r ] })
            .join(ch_blast_16s.result.map { meta, r -> [ meta.id, r ] })
            .map { _id, emeta, stats, kreport, blast -> [ emeta, stats, kreport, blast ] }
    )

    // Build ranked candidate species pool from all three pre-assembly tools
    ch_pool = build_candidates(
        ch_mash.distances
            .map  { meta, d -> [ meta.id, d ] }
            .join(ch_kraken.out.map    { meta, r -> [ meta.id, r ] })
            .join(ch_blast_16s.result.map { meta, r -> [ meta.id, r ] })
            .join(ch_meta_by_id)
            .map  { _id, distances, kreport, blast, emeta -> [ emeta, distances, kreport, blast ] }
    )

    // Download one RefSeq reference genome per candidate in the pool
    ch_refs = refseq_references(ch_pool.pool)

    // Multi-reference ANI confirmation with skani
    ch_skani = skani(
        ch_assembly_enriched
            .map  { meta, asm -> [ meta.id, asm ] }
            .join(ch_refs.references.map { meta, d -> [ meta.id, d ] })
            .join(ch_meta_by_id)
            .map  { _id, asm, refs_dir, emeta -> [ emeta, asm, refs_dir ] }
    )

    // BMGAP2 channels
    ch_bmgap2_amr = bmgap2_amr(
        ch_mlst.out.join(ch_pmga.out, by: 0)
    )
    ch_bmgap2_le  = bmgap2_locusextractor(ch_bmgap2_amr.out)
    ch_pre_report = bmgap2_bmscan(ch_bmgap2_le.out)

    // Species-specific analyses
    legsta(ch_assembly_enriched.filter          { meta, _f -> meta.mash_species == 'Legionella_pneumophila' })
    kleborate(ch_assembly_enriched.filter       { meta, _f -> meta.mash_genus   == 'Klebsiella' })
    shigatyper(ch_clean_enriched.filter         { meta, _f -> meta.mash_genus   == 'Shigella' })
    emm_typing(ch_clean_enriched.filter         { meta, _f -> meta.mash_species in ['Streptococcus_pyogenes', 'Streptococcus_dysgalactiae'] })
    seqsero2(ch_clean_enriched.filter           { meta, _f -> meta.mash_genus   == 'Salmonella' })
    serotypefinder(ch_clean_enriched.filter     { meta, _f -> meta.mash_species == 'Escherichia_coli' })
    plasmidfinder(ch_clean_enriched)
    seroba(ch_clean_enriched.filter             { meta, _f -> meta.mash_species == 'Streptococcus_pneumoniae' })
    pasty(ch_assembly_enriched.filter           { meta, _f -> meta.mash_species == 'Pseudomonas_aeruginosa' })
    kaptive_ab(ch_assembly_enriched.filter      { meta, _f -> meta.mash_species == 'Acinetobacter_baumannii' })
    kaptive_vp(ch_assembly_enriched.filter      { meta, _f -> meta.mash_species == 'Vibrio_parahaemolyticus' })
    lissero(ch_assembly_enriched.filter         { meta, _f -> meta.mash_species == 'Listeria_monocytogenes' })

    ch_optional_barrier =
        legsta.out.done
            .mix(kleborate.out.done, shigatyper.out.done, emm_typing.out.done,
                 seqsero2.out.done, serotypefinder.out.done, plasmidfinder.out.done,
                 seroba.out.done, pasty.out.done, kaptive_ab.out.done, kaptive_vp.out.done,
                 lissero.out.done,
                 ch_pre_report.map { meta, _f -> meta })
            .mix( ch_stats.map { meta, _s -> meta } )
            .mix( ch_skani.result.map { meta, _f -> meta } )
            .mix( ch_pool.pool.map { meta, _f -> meta } )
            .mix( ch_blast_16s.result.map { meta, _f -> meta } )
            .map { _id -> 1 }
            .collect()
            .map { _ids -> true }

    summary_report(
        ch_optional_barrier,
        ch_stats.map { _meta, stats -> stats }.collect(),
        ch_readssum.out.map         { _meta, rm   -> rm   }.collect(),
        ch_prokka.cds_txt.map       { _meta, ptxt -> ptxt }.collect(),
        ch_mlst.out.map             { _meta, mlst_file -> mlst_file }.collect(),
        ch_kraken_enriched.map      { _meta, kr   -> kr   }.collect(),
        ch_pmga.out.map             { _meta, pmga_file -> pmga_file }.collect().ifEmpty([]),
        ch_neisseria_txt,
        ch_hinfluenzae_txt,
        ch_aggregate.out.map        { _meta, f -> f }.collect().ifEmpty([]),
        ch_skani.result.map         { _meta, f -> f }.collect().ifEmpty([]),
        ch_blast_16s.result.map     { _meta, f -> f }.collect().ifEmpty([])
    )
}
