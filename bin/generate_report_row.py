#!/usr/bin/env python3
"""
generate_report_row.py — produce one headerless TSV data row per sample.

Output filename and column layout depend on mlst_scheme (read from pyoutputs):

  {sample_id}_row.tsv     — standard samples        (20 columns)
  {sample_id}_row_nm.tsv  — Neisseria meningitidis   (45 columns)
  {sample_id}_row_hi.tsv  — Haemophilus influenzae   (41 columns)

The header row is added once by summary_report.py.
"""

import argparse
import csv
import glob
import json
import os
import sys


# ---------------------------------------------------------------------------
# Sentinel helpers
# ---------------------------------------------------------------------------

NO_DATA = 'No data'


def normalize_le_value(val):
    if not val:
        return 'Not detected'
    if val in ['Not found', 'not found']:
        return 'Not detected'
    if val.startswith('Allele not identified') or val.startswith('Peptide not found'):
        return 'Not detected'
    return val


def find_gene(gene_dict, gene_name):
    for key, value in gene_dict.items():
        if value.get('Gene_name') == gene_name:
            return value
    for key, value in gene_dict.items():
        if key == gene_name or f"({gene_name})" in key or key.startswith(gene_name + " "):
            return value
    return None


# ---------------------------------------------------------------------------
# Species-specific serotype helpers  (non-meningitis)
# ---------------------------------------------------------------------------

def get_ecoli_serotype(mypath, sample_id):
    ecoli_json = os.path.join(mypath, "escherichia", "data.json")
    if not os.path.isfile(ecoli_json):
        return None
    try:
        with open(ecoli_json) as f:
            ecoli_data = json.load(f)
        sf_results = ecoli_data.get('serotypefinder', {}).get('results', {})
        o_types = list(dict.fromkeys(
            hit.get('serotype', '')
            for hit in sf_results.get('O_type', {}).values()
            if hit.get('serotype')
        ))
        h_types = list(dict.fromkeys(
            hit.get('serotype', '')
            for hit in sf_results.get('H_type', {}).values()
            if hit.get('serotype')
        ))
        o_str = '/'.join(sorted(set(o_types))) if o_types else 'NT'
        h_str = '/'.join(sorted(set(h_types))) if h_types else 'NT'
        return f"{o_str}:{h_str}"
    except Exception as e:
        print(f"Warning: Could not parse SerotypeFinder JSON for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_klebsiella_serotype(mypath, sample_id):
    kleb_tsv = os.path.join(mypath, "klebsiella", "kleborate-test-out.tsv")
    if not os.path.isfile(kleb_tsv):
        return None
    try:
        with open(kleb_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                st = row.get('ST', '').strip()
                return st if st and st != '-' else 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse Kleborate TSV for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_legionella_serotype(mypath, sample_id):
    legsta_txt = os.path.join(mypath, "legsta", "legsta_output.txt")
    if not os.path.isfile(legsta_txt):
        return None
    try:
        with open(legsta_txt) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                loci = [
                    row.get('SBT',   '-').strip(),
                    row.get('flaA',  '-').strip(),
                    row.get('pilE',  '-').strip(),
                    row.get('asd',   '-').strip(),
                    row.get('mip',   '-').strip(),
                    row.get('mompS', '-').strip(),
                    row.get('proA',  '-').strip(),
                    row.get('neuA',  '-').strip(),
                ]
                return ','.join(loci)
    except Exception as e:
        print(f"Warning: Could not parse Legsta output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_salmonella_serotype(mypath, sample_id):
    salm_dirs = glob.glob(os.path.join(mypath, "salmonella", "SeqSero_result_*"))
    if not salm_dirs:
        return None
    salm_tsv = os.path.join(salm_dirs[0], "SeqSero_result.tsv")
    if not os.path.isfile(salm_tsv):
        return None
    try:
        with open(salm_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                profile      = row.get('Predicted antigenic profile', '').strip()
                serotype_name = row.get('Predicted serotype', '').strip()
                if profile:
                    return f"{profile}({serotype_name})" if serotype_name else profile
                return 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse SeqSero2 TSV for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_gas_serotype(mypath, sample_id):
    gas_txt = os.path.join(mypath, "groupAstrep", "groupAstrep_result.txt")
    if not os.path.isfile(gas_txt):
        return None
    try:
        validated    = None
        nonvalidated = None
        with open(gas_txt) as f:
            for line in f:
                fields = line.strip().split()
                if len(fields) >= 3:
                    emm_raw = fields[2].replace('.sds', '')
                    if fields[1] == 'EMM_validated' and validated is None:
                        validated = emm_raw
                    elif fields[1] == 'EMM_nonValidated' and nonvalidated is None:
                        nonvalidated = emm_raw
        return validated or nonvalidated or 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse emm-typing-tool output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_shigella_serotype(mypath, sample_id):
    shigella_txt = os.path.join(mypath, "shigella", "shigella_output.txt")
    if not os.path.isfile(shigella_txt):
        return None
    try:
        with open(shigella_txt) as f:
            lines = f.readlines()
        found_header = False
        for line in lines:
            if found_header:
                fields = line.strip().split('\t')
                return fields[1].strip() if len(fields) >= 2 else 'Not detected'
            if line.startswith('sample\t'):
                found_header = True
        return 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse Shigatyper output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


# ---------------------------------------------------------------------------
# BMGAP2 data parsing
# ---------------------------------------------------------------------------

def parse_bmgap2_amr(mypath, sample_id, scheme, hinfluenzae_txt=None):
    bmgap2_data = {
        'penA_allele': NO_DATA,
        'penA_mutations': NO_DATA,
        'penA_phenotype': NO_DATA,
        'gyrA_allele': NO_DATA,
        'gyrA_mutations': NO_DATA,
        'gyrA_phenotype': NO_DATA,
        'parC_allele': NO_DATA,
        'parC_phenotype': NO_DATA,
        'rpoB_allele': NO_DATA,
        'rpoB_phenotype': NO_DATA,
        'ponA_allele': NO_DATA,
        'ponA_phenotype': NO_DATA,
        'predicted_resistance': NO_DATA,
        'bmgap2_species': NO_DATA,
        'bmgap2_mlst_st': NO_DATA,
        'bmgap2_mlst_cc': NO_DATA,
        'FHbp_variant': NO_DATA,
        'FHbp_subfamily': NO_DATA,
        'FHbp_peptide': NO_DATA,
        'NadA_variant': NO_DATA,
        'NhbA_peptide': NO_DATA,
        'vaccine_4CMenB_coverage': NO_DATA,
        'folA_allele': NO_DATA,
        'folA_phenotype': NO_DATA,
        'blaTEM1_status': NO_DATA,
        'blaROB1_status': NO_DATA,
    }

    amr_base = os.path.join(mypath, "bmgap2_amr")
    amr_json_files = glob.glob(os.path.join(amr_base, f"{sample_id}*amr_data.json"))

    if amr_json_files:
        try:
            with open(amr_json_files[0]) as f:
                amr_data = json.load(f)

            amr_genes = amr_data.get('amr_genes', {})

            pena = find_gene(amr_genes, 'penA')
            if pena:
                bmgap2_data['penA_allele'] = pena.get('allele', 'Not detected')
                known_muts = list(pena.get('known_mutations', {}).keys())
                bmgap2_data['penA_mutations'] = ';'.join(known_muts) if known_muts else 'None'
                if 'Penicillin' in amr_data.get('antimicrobics', {}).get('Penicillins', {}):
                    bmgap2_data['penA_phenotype'] = (
                        amr_data['antimicrobics']['Penicillins']['Penicillin']
                        .get('predicted_phenotype', 'Not detected')
                    )

            gyra = find_gene(amr_genes, 'gyrA')
            if gyra:
                bmgap2_data['gyrA_allele'] = gyra.get('allele', 'Not detected')
                known_muts = list(gyra.get('known_mutations', {}).keys())
                bmgap2_data['gyrA_mutations'] = ';'.join(known_muts) if known_muts else 'None'
                bmgap2_data['gyrA_phenotype'] = 'Susceptible' if not known_muts else 'Check'

            parc = find_gene(amr_genes, 'parC')
            if parc:
                bmgap2_data['parC_allele'] = parc.get('allele', 'Not detected')
                known_muts = list(parc.get('known_mutations', {}).keys())
                bmgap2_data['parC_phenotype'] = 'Susceptible' if not known_muts else 'Check'

            rpob = find_gene(amr_genes, 'rpoB')
            if rpob:
                bmgap2_data['rpoB_allele'] = rpob.get('allele', 'Not detected')
                known_muts = list(rpob.get('known_mutations', {}).keys())
                bmgap2_data['rpoB_phenotype'] = 'Susceptible' if not known_muts else 'Resistant'

            pona = find_gene(amr_genes, 'ponA')
            if pona:
                bmgap2_data['ponA_allele'] = pona.get('allele', 'Not detected')
                known_muts = list(pona.get('known_mutations', {}).keys())
                bmgap2_data['ponA_phenotype'] = 'Susceptible' if not known_muts else 'Check'

            # ftsI reuses penA_* slots for Hi (Hi has no penA)
            if scheme == 'hinfluenzae':
                ftsi = find_gene(amr_genes, 'ftsI')
                if ftsi:
                    bmgap2_data['penA_allele'] = ftsi.get('allele', 'Not detected')
                    known_muts = list(ftsi.get('known_mutations', {}).keys())
                    bmgap2_data['penA_mutations'] = ';'.join(known_muts) if known_muts else 'None'
                    if 'Penicillins' in amr_data.get('antimicrobics', {}):
                        for drug_data in amr_data['antimicrobics']['Penicillins'].values():
                            bmgap2_data['penA_phenotype'] = drug_data.get(
                                'predicted_phenotype', 'Not detected'
                            )
                            break
                    else:
                        bmgap2_data['penA_phenotype'] = (
                            'Susceptible' if not known_muts else 'Check'
                        )

            summary = amr_data.get('summary', {})
            bmgap2_data['predicted_resistance'] = summary.get('predicted_resistance', 'None')

            if scheme == 'hinfluenzae':
                fola = find_gene(amr_genes, 'folA')
                if fola:
                    bmgap2_data['folA_allele'] = fola.get('allele', 'Not detected')
                    known_muts = list(fola.get('known_mutations', {}).keys())
                    bmgap2_data['folA_phenotype'] = (
                        'Susceptible' if not known_muts else 'Resistant'
                    )

                blatem = find_gene(amr_genes, 'blaTEM-1')
                if blatem:
                    bmgap2_data['blaTEM1_status'] = blatem.get('status', 'Not detected')

                blarob = find_gene(amr_genes, 'blaROB-1')
                if blarob:
                    bmgap2_data['blaROB1_status'] = blarob.get('status', 'Not detected')

        except Exception as e:
            print(f"Warning: Could not parse AMR JSON for {sample_id}: {e}", file=sys.stderr)

    # LocusExtractor CSV
    le_base  = os.path.join(mypath, "bmgap2_locusextractor")
    le_dirs  = glob.glob(os.path.join(le_base, f"LE_*_{sample_id}_*"))
    le_csv_path = (
        os.path.join(le_dirs[0], "Results_text", f"molecular_data_{sample_id}.csv")
        if le_dirs else None
    )

    if le_csv_path and os.path.isfile(le_csv_path):
        try:
            with open(le_csv_path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if 'Assembly-1' not in row.get('Unique_ID', ''):
                        continue
                    if scheme == 'hinfluenzae':
                        bmgap2_data['bmgap2_mlst_st'] = (
                            row.get('Hi_MLST_ST', '') or NO_DATA
                        )
                        hi_st    = bmgap2_data['bmgap2_mlst_st']
                        hi_txt   = hinfluenzae_txt
                        bmgap2_data['bmgap2_mlst_cc'] = NO_DATA
                        if hi_txt and os.path.isfile(hi_txt) and hi_st not in [
                            NO_DATA, 'Not detected', 'New', 'NA', ''
                        ]:
                            try:
                                with open(hi_txt) as ht:
                                    for hrow in ht:
                                        hcols = hrow.rstrip().split('\t')
                                        if hcols[0] == hi_st:
                                            bmgap2_data['bmgap2_mlst_cc'] = (
                                                hcols[8] if len(hcols) >= 9 else 'Not detected'
                                            )
                                            break
                            except Exception:
                                pass
                    else:
                        bmgap2_data['bmgap2_mlst_st'] = (
                            row.get('Nm_MLST_ST', '') or NO_DATA
                        )
                        bmgap2_data['bmgap2_mlst_cc'] = (
                            row.get('Nm_MLST_cc', '') or NO_DATA
                        )

                    bmgap2_data['FHbp_variant']  = normalize_le_value(
                        row.get('FHbp_protein_subvariant_Novartis', '')
                    )
                    bmgap2_data['FHbp_subfamily'] = normalize_le_value(
                        row.get('FHbp_subfamily', '')
                    )
                    bmgap2_data['FHbp_peptide']  = normalize_le_value(
                        row.get('FHbp_protein_subvariant_Oxford', '')
                    )
                    bmgap2_data['NadA_variant']  = normalize_le_value(
                        row.get('NadA_Protein_subvariant_Novartis', '')
                    )
                    bmgap2_data['NhbA_peptide']  = normalize_le_value(
                        row.get('NhbA_Protein_subvariant_Novartis', '')
                    )

                    if scheme == 'hinfluenzae':
                        bmgap2_data['vaccine_4CMenB_coverage'] = 'Not applicable'
                    else:
                        has_fhbp = bmgap2_data['FHbp_variant'] not in [NO_DATA, 'Not detected']
                        has_nada = bmgap2_data['NadA_variant'] not in [NO_DATA, 'Not detected']
                        has_nhba = bmgap2_data['NhbA_peptide'] not in [NO_DATA, 'Not detected']
                        if has_fhbp and (has_nada or has_nhba):
                            bmgap2_data['vaccine_4CMenB_coverage'] = 'Likely'
                        elif has_fhbp:
                            bmgap2_data['vaccine_4CMenB_coverage'] = 'Possible'
                        else:
                            bmgap2_data['vaccine_4CMenB_coverage'] = 'Unlikely'
                    break

        except Exception as e:
            print(
                f"Warning: Could not parse LocusExtractor CSV for {sample_id}: {e}",
                file=sys.stderr
            )

    # BMScan JSON
    bmscan_base = os.path.join(mypath, "bmgap2_bmscan")
    if os.path.isdir(bmscan_base):
        bmscan_json_files = glob.glob(os.path.join(bmscan_base, "species_analysis_*.json"))
        if bmscan_json_files:
            try:
                with open(bmscan_json_files[0]) as f:
                    bmscan_data = json.load(f)
                for sample_data in bmscan_data.values():
                    if 'mash_results' in sample_data:
                        bmgap2_data['bmgap2_species'] = (
                            sample_data['mash_results'].get('species', '-')
                        )
                        break
            except Exception as e:
                print(
                    f"Warning: Could not parse BMScan JSON for {sample_id}: {e}",
                    file=sys.stderr
                )

    return bmgap2_data


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate one headerless TSV row for a single sample."
    )
    parser.add_argument("pyoutputs",   help="{sample_id}_pyoutputs.txt")
    parser.add_argument("sample_id",   help="Sample ID")
    parser.add_argument("mypath",      help="Full path to sample output directory")
    parser.add_argument(
        "--hinfluenzae_txt",
        default=None,
        help="Path to hinfluenzae MLST CC lookup table (required for Hi samples)"
    )
    parser.add_argument(
        "--meningitis",
        choices=["true", "false"],
        default="false",
        help="Enable BMGAP2 meningitis data parsing"
    )
    args = parser.parse_args()

    sample_id  = args.sample_id
    mypath     = args.mypath
    meningitis = args.meningitis == "true"

    # --- Read pyoutputs CSV ---
    with open(args.pyoutputs) as f:
        for line in f:
            cells = line.rstrip().split(",")

    scheme = cells[16]

    # --- Serotype notes (PMGA sta file, meningitis species only) ---
    serotype_notes = '-'
    pmga_sta_file  = os.path.join(mypath, "pmga", f"{sample_id}sta.txt")
    if os.path.isfile(pmga_sta_file):
        try:
            with open(pmga_sta_file) as f:
                lines = f.readlines()
            if len(lines) > 1:
                data_line = lines[1].rstrip().split('\t')
                if len(data_line) >= 5:
                    serotype_notes = data_line[4]
        except Exception as e:
            print(
                f"Warning: Could not parse PMGA sta file for {sample_id}: {e}",
                file=sys.stderr
            )

    # --- Serotype ---
    serotype_value = cells[22]  # PMGA prediction — nm_serogroup / hi_serotype for Nm/Hi
    if scheme not in ('neisseria', 'hinfluenzae'):
        # Species-specific typing tools override PMGA for standard samples
        candidates = [
            get_ecoli_serotype(mypath, sample_id),
            get_klebsiella_serotype(mypath, sample_id),
            get_legionella_serotype(mypath, sample_id),
            get_salmonella_serotype(mypath, sample_id),
            get_gas_serotype(mypath, sample_id),
            get_shigella_serotype(mypath, sample_id),
        ]
        for c in candidates:
            if c is not None:
                serotype_value = c
                break
        serotype_notes = NO_DATA

    # --- BMGAP2 data (Nm/Hi only, requires meningitis=true) ---
    bmgap2_data = None
    if meningitis and scheme in ('neisseria', 'hinfluenzae'):
        bmgap2_data = parse_bmgap2_amr(mypath, sample_id, scheme, hinfluenzae_txt=args.hinfluenzae_txt)

    def bd(key):
        return bmgap2_data[key] if bmgap2_data else NO_DATA

    # --- Build species-appropriate row ---
    pmga_species = cells[21] if cells[21] else NO_DATA

    common = [
        sample_id,
        cells[13], cells[11], cells[12], cells[14],
        cells[5], cells[6], cells[7], cells[8], cells[9], cells[10],
        cells[15],
        cells[0] + '_' + cells[1], cells[3], cells[2],
        cells[20], cells[19],
    ]

    if scheme == 'neisseria':
        row = common + [
            cells[16], cells[17], cells[18],             # mlst_scheme, mlst_st, mlst_cc
            pmga_species,                                # pmga_species
            serotype_value,                              # nm_serogroup
            serotype_notes,
            bd('bmgap2_species'), bd('bmgap2_mlst_st'), bd('bmgap2_mlst_cc'), bd('predicted_resistance'),
            bd('penA_allele'), bd('penA_mutations'), bd('penA_phenotype'),
            bd('gyrA_allele'), bd('gyrA_mutations'), bd('gyrA_phenotype'),
            bd('parC_allele'), bd('parC_phenotype'),
            bd('rpoB_allele'), bd('rpoB_phenotype'),
            bd('ponA_allele'), bd('ponA_phenotype'),
            bd('FHbp_variant'), bd('FHbp_subfamily'), bd('FHbp_peptide'),
            bd('NadA_variant'), bd('NhbA_peptide'),
            bd('vaccine_4CMenB_coverage'),
        ]
        out_file = f"{sample_id}_row_nm.tsv"

    elif scheme == 'hinfluenzae':
        row = common + [
            cells[16], cells[17], cells[18],             # mlst_scheme, mlst_st, mlst_cc
            pmga_species,                                # pmga_species
            serotype_value,                              # hi_serotype
            serotype_notes,
            bd('bmgap2_species'), bd('bmgap2_mlst_st'), bd('bmgap2_mlst_cc'), bd('predicted_resistance'),
            bd('penA_allele'), bd('penA_mutations'), bd('penA_phenotype'),  # ftsI data in penA slots
            bd('gyrA_allele'), bd('gyrA_mutations'), bd('gyrA_phenotype'),
            bd('parC_allele'), bd('parC_phenotype'),
            bd('rpoB_allele'), bd('rpoB_phenotype'),
            bd('folA_allele'), bd('folA_phenotype'),
            bd('blaTEM1_status'), bd('blaROB1_status'),
        ]
        out_file = f"{sample_id}_row_hi.tsv"

    else:
        row = common + [
            cells[16], cells[17],                        # mlst_scheme, mlst_st
            serotype_value,
        ]
        out_file = f"{sample_id}_row.tsv"

    with open(out_file, 'w') as fh:
        fh.write('\t'.join(str(v) for v in row) + '\n')

    print(f"Wrote {out_file}")


if __name__ == "__main__":
    main()
