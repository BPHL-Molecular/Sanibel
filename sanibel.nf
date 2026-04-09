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
include { bbtools }               from './modules/bbtools.nf'
include { fastqc2 }               from './modules/fastqc2.nf'
include { multiqc }               from './modules/multiqc.nf'
include { mash }                  from './modules/mash.nf'
include { unicycler }             from './modules/unicycler.nf'
include { quast }                 from './modules/quast.nf'
include { kraken }                from './modules/kraken.nf'
include { parse_assembly }        from './modules/parse_assembly.nf'
include { readssum }              from './modules/readssum.nf'
include { parse_reads }           from './modules/parse_reads.nf'
include { prokka }                from './modules/prokka.nf'
include { amrfinder }             from './modules/amrfinder.nf'
include { mlst }                  from './modules/mlst.nf'
include { pmga }                  from './modules/pmga.nf'
include { parse_typing }          from './modules/parse_typing.nf'
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
include { generate_row }          from './modules/generate_row.nf'
include { summary_report }        from './modules/summary_report.nf'

workflow {
    log.info """
    Sanibel — Bacterial WGS Analysis Pipeline
    ==========================================================================
    input dir   : ${params.input}
    output dir  : ${params.output}
    bmgap2 db   : ${params.bmgap2_db}
    ==========================================================================
    """

    // Input channel
    ch_reads = channel.fromFilePairs(
        ["${params.input}/*_{1,2}.fastq.gz",
         "${params.input}/*_R{1,2}_*.fastq.gz"],
        checkIfExists: true
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
          "https://raw.githubusercontent.com/tseemann/mlst/master/db/pubmlst/haemophilus_influenzae/haemophilus_influenzae.txt" ]
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
    ch_clean   = bbtools(ch_trimmed.reads)
    ch_fastqc2 = fastqc2(ch_clean.reads)

    // Per-sample MultiQC combining raw and clean FastQC reports
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

    // Assembly statistics parsing
    ch_parse_assembly = parse_assembly(
        ch_mash.distances
            .join(ch_quast.report, by: 0)
    )

    // Read metrics
    ch_readssum = readssum(
        ch_parse_assembly.out
            .join(ch_clean.reads, by: 0)
    )
    ch_parse_reads = parse_reads(ch_readssum.out)

    // Annotation and typing
    ch_prokka = prokka(
        ch_assembly.assembly
            .join(ch_parse_reads.out, by: 0)
    )
    amrfinder(ch_assembly.assembly)
    ch_mlst = mlst(ch_assembly.assembly)
    ch_pmga = pmga(ch_assembly.assembly.join(ch_mlst.out, by: 0))

    // Typing data collection
    ch_parse_typing = parse_typing(
        ch_prokka.cds
            .join(ch_mlst.out, by: 0)
            .join(ch_kraken.out, by: 0)
            .join(ch_pmga.out, by: 0)
            .combine(ch_neisseria_txt)
            .combine(ch_hinfluenzae_txt)
    )

    // BMGAP2 modules — Python scripts self-gate on Nm/Hi species
    ch_bmgap2_amr = bmgap2_amr(ch_parse_typing.out)
    ch_bmgap2_le  = bmgap2_locusextractor(ch_bmgap2_amr.out)
    ch_pre_report = bmgap2_bmscan(ch_bmgap2_le.out)

    // Species-specific analyses
    ch_pre_assembly = ch_pre_report.join(ch_assembly.assembly, by: 0)
    ch_pre_reads    = ch_pre_report.join(ch_clean.reads,       by: 0)

    legsta(ch_pre_assembly)
    kleborate(ch_pre_assembly)
    shigatyper(ch_pre_reads)
    emm_typing(ch_pre_reads)
    seqsero2(ch_pre_reads)
    serotypefinder(ch_pre_reads)
    plasmidfinder(ch_pre_reads)
    seroba(ch_pre_reads)
    pasty(ch_pre_assembly)
    kaptive_ab(ch_pre_assembly)
    kaptive_vp(ch_pre_assembly)

    // Per-sample report row
    ch_species_done =
        legsta.done
            .mix(shigatyper.done, emm_typing.done, kleborate.done,
                 seqsero2.done, serotypefinder.done, plasmidfinder.done,
                 seroba.done, pasty.done, kaptive_ab.done, kaptive_vp.done)
            .map { meta -> [meta.id, 1] }
            .groupTuple(size: 11)

    ch_generate_row = generate_row(
        ch_pre_report
            .map { meta, pyoutputs -> [meta.id, meta, pyoutputs] }
            .join(ch_species_done)
            .map { _id, meta, pyoutputs, _done -> [meta, pyoutputs] }
            .combine(ch_hinfluenzae_txt)
    )

    summary_report(
        ch_generate_row.row
            .map { _meta, row -> row }
            .flatten()
            .collect()
    )
}
