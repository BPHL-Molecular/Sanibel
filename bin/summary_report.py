#!/usr/bin/env python3
"""
summary_report.py — collect per-sample row files and write three summary reports.

Output files (only created if samples of that type were processed):
    sum_report.txt     — standard samples        (20 columns)
    sum_report_nm.txt  — Neisseria meningitidis   (45 columns)
    sum_report_hi.txt  — Haemophilus influenzae   (41 columns)
"""

import glob
import os
import sys


HEADER_STANDARD = [
    'sampleID',
    'num_clean_reads',
    'avg_readlength',
    'avg_read_qual',
    'est_coverage',
    'num_contigs',
    'longest_contig',
    'N50',
    'L50',
    'total_length',
    'gc_content',
    'annotated_cds',
    'speciesID_mash',
    'nearest_neighbor_mash',
    'mash_distance',
    'speciesID_kraken',
    'kraken_percent',
    'mlst_scheme',
    'mlst_st',
    'serotype',
]

HEADER_NM = [
    'sampleID',
    'num_clean_reads',
    'avg_readlength',
    'avg_read_qual',
    'est_coverage',
    'num_contigs',
    'longest_contig',
    'N50',
    'L50',
    'total_length',
    'gc_content',
    'annotated_cds',
    'speciesID_mash',
    'nearest_neighbor_mash',
    'mash_distance',
    'speciesID_kraken',
    'kraken_percent',
    'mlst_scheme',
    'mlst_st',
    'mlst_cc',
    'pmga_species',
    'nm_serogroup',
    'serotype_notes',
    'bmgap2_species',
    'bmgap2_mlst_st',
    'bmgap2_mlst_cc',
    'predicted_resistance',
    'penA_allele',
    'penA_mutations',
    'penA_phenotype',
    'gyrA_allele',
    'gyrA_mutations',
    'gyrA_phenotype',
    'parC_allele',
    'parC_phenotype',
    'rpoB_allele',
    'rpoB_phenotype',
    'ponA_allele',
    'ponA_phenotype',
    'FHbp_variant',
    'FHbp_subfamily',
    'FHbp_peptide',
    'NadA_variant',
    'NhbA_peptide',
    'vaccine_4CMenB_coverage',
]

HEADER_HI = [
    'sampleID',
    'num_clean_reads',
    'avg_readlength',
    'avg_read_qual',
    'est_coverage',
    'num_contigs',
    'longest_contig',
    'N50',
    'L50',
    'total_length',
    'gc_content',
    'annotated_cds',
    'speciesID_mash',
    'nearest_neighbor_mash',
    'mash_distance',
    'speciesID_kraken',
    'kraken_percent',
    'mlst_scheme',
    'mlst_st',
    'mlst_cc',
    'pmga_species',
    'hi_serotype',
    'serotype_notes',
    'bmgap2_species',
    'bmgap2_mlst_st',
    'bmgap2_mlst_cc',
    'predicted_resistance',
    'ftsI_allele',
    'ftsI_mutations',
    'ftsI_phenotype',
    'gyrA_allele',
    'gyrA_mutations',
    'gyrA_phenotype',
    'parC_allele',
    'parC_phenotype',
    'rpoB_allele',
    'rpoB_phenotype',
    'folA_allele',
    'folA_phenotype',
    'blaTEM1_status',
    'blaROB1_status',
]


def write_report(outfile, header, row_files, suffix):
    rows = {}
    for f in row_files:
        sample_id = f[:-len(suffix)]
        rows[sample_id] = f

    lines = ['\t'.join(header)]
    for sample_id in sorted(rows.keys()):
        with open(rows[sample_id]) as fh:
            lines.append(fh.read().strip())

    with open(outfile, 'w') as fh:
        fh.write(os.linesep.join(lines) + os.linesep)

    print(f"summary_report.py: wrote {outfile} ({len(rows)} sample(s))")


def main():
    standard_files = glob.glob("*_row.tsv")
    nm_files       = glob.glob("*_row_nm.tsv")
    hi_files       = glob.glob("*_row_hi.tsv")

    if not (standard_files or nm_files or hi_files):
        print("summary_report.py: no row files found in working directory.", file=sys.stderr)
        sys.exit(1)

    if standard_files:
        write_report("sum_report.txt",    HEADER_STANDARD, standard_files, "_row.tsv")
    if nm_files:
        write_report("sum_report_nm.txt", HEADER_NM,       nm_files,       "_row_nm.tsv")
    if hi_files:
        write_report("sum_report_hi.txt", HEADER_HI,       hi_files,       "_row_hi.tsv")


if __name__ == "__main__":
    main()
