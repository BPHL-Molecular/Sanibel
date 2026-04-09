#!/usr/bin/env python3
"""
collect_sample_data.py — unified pyoutputs accumulation script.

Subcommands (must be run in order):
  parse-assembly  <mash_top10_tab> <quast_report_tsv> <sample_id>_pyoutputs.txt
  parse-reads     <read_metrics_txt> <sample_id>_pyoutputs.txt
  parse-typing    <mypath> <sample_id> <sample_id>_pyoutputs.txt
"""

import argparse
import os
import re
import sys

import pandas as pd


# ---------------------------------------------------------------------------
# parse-assembly  (replaces pyTask1 logic)
# ---------------------------------------------------------------------------

def cmd_parse_assembly(args):
    mash_tab    = args.mash_tab
    quast_tsv   = args.quast_tsv
    output_file = args.output_file

    with open(mash_tab, 'r') as mash:
        top_hit = mash.readline()
        top_hit = str(top_hit)
        gn      = re.sub('.*-\.-', '', top_hit)
        cells   = gn.split()
        acell   = cells[0]
        acell   = acell.lstrip("_")
        agn     = re.split(r'^([^_]*_[^_]*)(_|\.).*$', acell)[1]
        genus   = agn.split('_')[0]
        species = agn.split('_')[1]
        distance   = top_hit.split()[2]
        accession  = top_hit.split("-")[5]

    df         = pd.read_table(quast_tsv, sep="\t")
    assem      = list(df.columns)[1]
    contigs    = df[assem][12].astype(int)
    long_contig = df[assem][13].astype(int)
    n50        = df[assem][16].astype(int)
    l50        = df[assem][18].astype(int)
    genome     = df[assem][14].astype(int)
    gc         = df[assem][15].astype(int)

    with open(output_file, "w") as f:
        f.write(
            f"{genus},{species},{distance},{accession},{assem},"
            f"{contigs},{long_contig},{n50},{l50},{genome},{gc}"
        )


# ---------------------------------------------------------------------------
# parse-reads  (replaces pyTask2 logic)
# ---------------------------------------------------------------------------

def cmd_parse_reads(args):
    metrics_file = args.metrics_file
    output_file  = args.output_file

    with open(metrics_file, 'r') as metrics:
        firstline  = metrics.readline().rstrip().split()
        secondline = metrics.readline().rstrip().split()
        rm = dict(zip(firstline, secondline))
        avg_read_len = rm['avgReadLength']
        avg_qual     = rm['avgQuality']
        num_read     = rm['numReads']
        cov          = rm['coverage']

    with open(output_file, "a") as f:
        f.write(f",{avg_read_len},{avg_qual},{num_read},{cov}")


# ---------------------------------------------------------------------------
# parse-typing  (replaces pyTask3 logic)
# ---------------------------------------------------------------------------

def cmd_parse_typing(args):
    mypath      = args.mypath
    sample_id   = args.sample_id
    output_file = args.output_file

    filepath1  = args.prokka_txt      # path to prokka {sample}.txt (CDS count)
    filepath2  = args.mlst_file       # path to {sample}.mlst
    filepath3  = args.kraken_report   # path to Kraken2 .report
    filepath5a = args.neisseria_txt   # path to neisseria MLST CC table
    filepath5b = args.hinfluenzae_txt # path to hinfluenzae MLST CC table

    cds = ''
    with open(filepath1, 'r') as genes:
        for line in genes:
            line    = line.rstrip()
            content = line.split()
            if content and content[0] == 'CDS:':
                cds = content[1]

    scheme = ''
    st     = ''
    cc     = ''
    with open(filepath2, 'r') as mlst:
        for line in mlst:
            out    = line.rstrip().split()
            scheme = out[1]
            st     = out[2]
            cc     = ""

            if scheme == "neisseria":
                with open(filepath5a, 'r') as mlst_table:
                    for row in mlst_table:
                        cols = row.rstrip().split("\t")
                        if st == cols[0]:
                            cc = cols[8] if len(cols) == 9 else "NA"
                            break

            if scheme == "hinfluenzae":
                with open(filepath5b, 'r') as mlst_table:
                    for row in mlst_table:
                        cols = row.rstrip().split("\t")
                        if st == cols[0]:
                            cc = cols[8] if len(cols) == 9 else "NA"
                            break

    percent = ''
    tax     = ''
    with open(filepath3, 'r') as kreport:
        for l in kreport:
            l_parse  = l.lstrip().rstrip().split("\t")
            percent  = l_parse[0]
            tax_level = l_parse[3]
            tax      = l_parse[5].lstrip()
            if tax_level == 'S':
                break

    pgspecies = ''
    pgpredict = ''
    if os.path.getsize(args.pmga_txt) > 0:
        with open(args.pmga_txt, 'r') as pgmalines:
            pglines  = pgmalines.readlines()
            pgcells  = pglines[1].strip().split("\t")
            pgpredict = pgcells[2]
            pgspecies = pgcells[1]

    with open(output_file, "a") as f:
        f.write(
            f",{cds},{scheme},{st},{cc},{percent},{tax},{pgspecies},{pgpredict}"
        )


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Accumulate per-sample data into a pyoutputs CSV file."
    )
    subparsers = parser.add_subparsers(dest="subcommand", required=True)

    # parse-assembly
    p_assem = subparsers.add_parser(
        "parse-assembly",
        help="Create pyoutputs file from mash + quast outputs (pyTask1 replacement)"
    )
    p_assem.add_argument("mash_tab",    help="Path to {sample}_distances_top10.tab")
    p_assem.add_argument("quast_tsv",   help="Path to quast_results/report.tsv")
    p_assem.add_argument("output_file", help="Output file to create (e.g. {sample}_pyoutputs.txt)")
    p_assem.set_defaults(func=cmd_parse_assembly)

    # parse-reads
    p_reads = subparsers.add_parser(
        "parse-reads",
        help="Append read metrics to pyoutputs file (pyTask2 replacement)"
    )
    p_reads.add_argument("metrics_file", help="Path to {sample}_readMetrics.txt")
    p_reads.add_argument("output_file",  help="Existing pyoutputs file to append to")
    p_reads.set_defaults(func=cmd_parse_reads)

    # parse-typing
    p_typing = subparsers.add_parser(
        "parse-typing",
        help="Append MLST/Kraken/PMGA results to pyoutputs file (pyTask3 replacement)"
    )
    p_typing.add_argument("mypath",           help="Full path to sample output directory (used for PMGA invocation)")
    p_typing.add_argument("sample_id",        help="Sample ID")
    p_typing.add_argument("output_file",      help="Existing pyoutputs file to append to")
    p_typing.add_argument("--prokka_txt",     required=True, help="Path to prokka {sample}.txt file (CDS count)")
    p_typing.add_argument("--mlst_file",      required=True, help="Path to {sample}.mlst file")
    p_typing.add_argument("--kraken_report",  required=True, help="Path to Kraken2 report file")
    p_typing.add_argument("--pmga_txt",        required=True, help="Path to PMGA {sample}sta.txt output (empty file if not Neisseria/H.influenzae)")
    p_typing.add_argument("--neisseria_txt",   required=True, help="Path to neisseria MLST CC lookup table")
    p_typing.add_argument("--hinfluenzae_txt", required=True, help="Path to hinfluenzae MLST CC lookup table")
    p_typing.set_defaults(func=cmd_parse_typing)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
