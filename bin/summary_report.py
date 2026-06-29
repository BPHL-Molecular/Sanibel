#!/usr/bin/env python3
"""
summary_report.py — Build Sanibel summary reports from individual tool outputs.

"""

import argparse
import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sanibel_taxonomy import (
    SYNONYMOUS_PAIRS, genus_of, detect_contamination, extract_contam_candidates,
)


NO_DATA = 'No data'

def normalize_le_value(val):
    if not val:
        return 'Not detected'
    if val in ['Not found', 'not found']:
        return 'Not detected'
    if val.startswith('Allele not identified') or val.startswith('Peptide not found'):
        return 'Not detected'
    if val.startswith('Incomplete ORF') or val.startswith('incomplete ORF'):
        return 'Not detected'
    if val.startswith('New-BLASTonly') or val.startswith('New-PCR'):
        return 'New allele'
    return val


def find_gene(gene_dict, gene_name):
    for key, value in gene_dict.items():
        if value.get('Gene_name') == gene_name:
            return value
    for key, value in gene_dict.items():
        if key == gene_name or f"({gene_name})" in key or key.startswith(gene_name + ' '):
            return value
    return None


def parse_assembly_stats(filepath):
    with open(filepath) as f:
        fields = f.read().strip().split(',')
    return {
        'genus':          fields[0],
        'species':        fields[1],
        'mash_distance':  fields[2],
        'accession':      fields[3],
        'num_contigs':    fields[5],
        'longest_contig': fields[6],
        'n50':            fields[7],
        'l50':            fields[8],
        'total_length':   fields[9],
        'gc_content':     fields[10],
    }


def parse_read_metrics(filepath):
    with open(filepath) as f:
        header = f.readline().strip().split()
        values = f.readline().strip().split()
    rm = dict(zip(header, values))
    return {
        'avg_read_len': rm.get('avgReadLength', NO_DATA),
        'avg_qual':     rm.get('avgQuality',    NO_DATA),
        'num_reads':    rm.get('numReads',      NO_DATA),
        'coverage':     rm.get('coverage',      NO_DATA),
    }


def parse_prokka_txt(filepath):
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if parts and parts[0] == 'CDS:':
                return parts[1]
    return NO_DATA


def parse_mlst(filepath, neisseria_txt=None, hinfluenzae_txt=None):
    scheme = st = cc = ''
    with open(filepath) as f:
        for line in f:
            out = line.strip().split()
            if len(out) >= 3:
                scheme = out[1]
                st     = out[2]
                cc     = ''
                if scheme in ('-', '') :
                    scheme = 'Not detected'
                if st in ('-', ''):
                    st = 'Not detected'
                if scheme == 'neisseria' and neisseria_txt and os.path.isfile(neisseria_txt):
                    try:
                        with open(neisseria_txt) as tbl:
                            for row in tbl:
                                cols = row.strip().split('\t')
                                if cols[0] == st:
                                    cc = cols[8] if len(cols) >= 9 else 'NA'
                                    break
                    except Exception:
                        pass
                elif scheme == 'hinfluenzae' and hinfluenzae_txt and os.path.isfile(hinfluenzae_txt):
                    try:
                        with open(hinfluenzae_txt) as tbl:
                            for row in tbl:
                                cols = row.strip().split('\t')
                                if cols[0] == st:
                                    cc = cols[8] if len(cols) >= 9 else 'NA'
                                    break
                    except Exception:
                        pass
            break
    return {'scheme': scheme, 'st': st, 'cc': cc}


def parse_kraken_report(filepath):
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6 and parts[3] == 'S':
                return {'percent': parts[0].strip(), 'species': parts[5].strip()}
    return {'percent': NO_DATA, 'species': NO_DATA}


def parse_pmga(filepath):
    result = {'species': '', 'prediction': '', 'serotype_notes': '-'}
    try:
        if os.path.getsize(filepath) > 0:
            with open(filepath) as f:
                lines = f.readlines()
            if len(lines) > 1:
                data = lines[1].strip().split('\t')
                result['species']        = data[1] if len(data) > 1 else ''
                result['prediction']     = data[2] if len(data) > 2 else ''
                result['serotype_notes'] = data[4] if len(data) > 4 else '-'
    except Exception as e:
        print(f"Warning: Could not parse PMGA file {filepath}: {e}", file=sys.stderr)
    return result


def parse_aggregate_species(filepath):
    empty = {'genus': NO_DATA, 'contamination_flag': NO_DATA}
    try:
        with open(filepath) as f:
            f.readline()  # skip header
            line = f.readline().strip()
        if not line:
            return empty
        parts = line.split('\t')
        genus       = parts[0] if len(parts) > 0 else NO_DATA
        contam_flag = parts[4] if len(parts) > 4 else NO_DATA
        return {'genus': genus, 'contamination_flag': contam_flag}
    except Exception:
        return empty


def parse_skani(filepath):
    EMPTY = {'ani': 'ANI < 80%', 'confirmed_species': 'NO ID', 'align_fraction': NO_DATA, 'reference': NO_DATA}
    try:
        rows = []
        with open(filepath) as f:
            f.readline()
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 5:
                    continue
                try:
                    ani = float(parts[2])
                except ValueError:
                    continue
                rows.append((ani, parts))
        if not rows:
            return EMPTY
        rows.sort(key=lambda x: x[0], reverse=True)
        best_ani, best_parts = rows[0]
        ref_basename = os.path.basename(best_parts[0])
        if ref_basename.endswith('.fna'):
            ref_basename = ref_basename[:-4]
        # References are named <Genus_species>__<ACCESSION>.fna (multiple strains per
        # candidate); recover the clean species label and accession from the filename,
        # falling back to the legacy single-file naming + Ref_name column otherwise.
        if '__' in ref_basename:
            species_part, acc_part = ref_basename.split('__', 1)
            confirmed_species = species_part if species_part else NO_DATA
            skani_ref_acc     = acc_part if acc_part else NO_DATA
        else:
            confirmed_species = ref_basename if ref_basename else NO_DATA
            ref_name_col      = best_parts[5].strip() if len(best_parts) > 5 else ''
            skani_ref_acc     = ref_name_col.split()[0] if ref_name_col else NO_DATA
        align_fraction = best_parts[4] if len(best_parts) > 4 else NO_DATA
        return {
            'ani':               f"{best_ani:.3f}",
            'confirmed_species': confirmed_species,
            'align_fraction':    align_fraction,
            'reference':         skani_ref_acc if skani_ref_acc else NO_DATA,
        }
    except Exception:
        return EMPTY


def parse_blast16s_result(filepath, anchor_genera=None):
    MIN_LENGTH = 400
    MIN_PIDENT = 97.0
    qualifying = []
    try:
        with open(filepath) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 13 or parts[0] == 'qseqid':
                    continue
                try:
                    pident   = float(parts[2])
                    length   = int(parts[3])
                    bitscore = float(parts[11])
                except ValueError:
                    continue
                if length < MIN_LENGTH or pident < MIN_PIDENT:
                    continue
                stitle = parts[12].strip()
                genus  = stitle.split()[0] if stitle else ''
                qualifying.append((pident, length, genus, bitscore))
    except Exception:
        pass
    if not qualifying:
        return {'pident': NO_DATA, 'tophit': NO_DATA}

    # Top hit = highest pident; closest to ~1500 bp as tiebreaker; highest bitscore.
    qualifying.sort(key=lambda r: (-r[0], abs(r[1] - 1500), -r[3]))

    # Prefer the highest-identity hit whose genus agrees with Mash or Kraken;
    # fall back to the overall top hit if 16S agrees with neither.
    anchor = {g.lower() for g in (anchor_genera or []) if g and g not in (NO_DATA, 'Unknown')}
    chosen = next((r for r in qualifying if r[2].lower() in anchor), qualifying[0]) if anchor \
             else qualifying[0]

    pident, _len, genus, _bs = chosen
    return {'pident': f"{pident:.3f}", 'tophit': f"{genus} spp." if genus else NO_DATA}


# Per-file parsers — species-specific

def get_ecoli_serotype(sample_dir, sample_id):
    ecoli_json = os.path.join(sample_dir, 'serotypefinder', 'data.json')
    if not os.path.isfile(ecoli_json):
        return None
    try:
        with open(ecoli_json) as f:
            data = json.load(f)
        sf = data.get('serotypefinder', {}).get('results', {})
        o_types = list(dict.fromkeys(
            hit.get('serotype', '') for hit in sf.get('O_type', {}).values() if hit.get('serotype')
        ))
        h_types = list(dict.fromkeys(
            hit.get('serotype', '') for hit in sf.get('H_type', {}).values() if hit.get('serotype')
        ))
        o_str = '/'.join(sorted(set(o_types))) if o_types else 'NT'
        h_str = '/'.join(sorted(set(h_types))) if h_types else 'NT'
        return f"{o_str}:{h_str}"
    except Exception as e:
        print(f"Warning: Could not parse SerotypeFinder JSON for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_klebsiella_serotype(sample_dir, sample_id):
    kleb_tsv = os.path.join(sample_dir, 'kleborate', 'kleborate_out',
                            'klebsiella_pneumo_complex_output.txt')
    if not os.path.isfile(kleb_tsv):
        return None
    try:
        with open(kleb_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                st = row.get('ST', '').strip()
                return st if st and st != '-' else 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse Kleborate output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_legionella_serotype(sample_dir, sample_id):
    legsta_txt = os.path.join(sample_dir, 'legsta', 'legsta_output.txt')
    if not os.path.isfile(legsta_txt):
        return None
    try:
        with open(legsta_txt) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                loci = [
                    row.get('SBT',   '-').strip(), row.get('flaA',  '-').strip(),
                    row.get('pilE',  '-').strip(), row.get('asd',   '-').strip(),
                    row.get('mip',   '-').strip(), row.get('mompS', '-').strip(),
                    row.get('proA',  '-').strip(), row.get('neuA',  '-').strip(),
                ]
                return ','.join(loci)
    except Exception as e:
        print(f"Warning: Could not parse Legsta output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_salmonella_serotype(sample_dir, sample_id):
    salm_dirs = glob.glob(os.path.join(sample_dir, 'seqsero2', 'SeqSero_result_*'))
    if not salm_dirs:
        return None
    salm_tsv = os.path.join(salm_dirs[0], 'SeqSero_result.tsv')
    if not os.path.isfile(salm_tsv):
        return None
    try:
        with open(salm_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                profile  = row.get('Predicted antigenic profile', '').strip()
                serotype = row.get('Predicted serotype', '').strip()
                if profile:
                    return f"{profile}({serotype})" if serotype else profile
                return 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse SeqSero2 TSV for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_gas_serotype(sample_dir, sample_id):
    gas_txt = os.path.join(sample_dir, 'emm_typing', 'groupAstrep_result.txt')
    if not os.path.isfile(gas_txt):
        return None
    try:
        validated = nonvalidated = None
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


def get_shigella_serotype(sample_dir, sample_id):
    sh_txt = os.path.join(sample_dir, 'shigatyper', 'shigatyper_output.txt')
    if not os.path.isfile(sh_txt):
        return None
    try:
        found_header = False
        with open(sh_txt) as f:
            for line in f:
                if found_header:
                    fields = line.strip().split('\t')
                    return fields[1].strip() if len(fields) >= 2 else 'Not detected'
                if line.startswith('sample\t'):
                    found_header = True
        return 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse Shigatyper output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_pneumococcal_serotype(sample_dir, sample_id):
    pred_csv = os.path.join(sample_dir, 'seroba', sample_id, 'pred.csv')
    if not os.path.isfile(pred_csv):
        return None
    try:
        with open(pred_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                serotype = row.get('Serotype', '').strip()
                return serotype if serotype else 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse Seroba pred.csv for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def _parse_kaptive_txt(filepath):
    """Return (locus, type) from a kaptive v3 TSV output file."""
    try:
        with open(filepath) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                locus = row.get('Best match locus', '').strip()
                ktype = row.get('Best match type', '').strip()
                return locus, ktype
    except Exception:
        pass
    return '', ''


def get_acinetobacter_serotype(sample_dir, sample_id):
    k_txt  = os.path.join(sample_dir, 'kaptive_ab', f'{sample_id}_ab_k.txt')
    oc_txt = os.path.join(sample_dir, 'kaptive_ab', f'{sample_id}_ab_oc.txt')
    if not os.path.isfile(k_txt) and not os.path.isfile(oc_txt):
        return None
    try:
        k_locus,  k_type  = _parse_kaptive_txt(k_txt)  if os.path.isfile(k_txt)  else ('', '')
        oc_locus, oc_type = _parse_kaptive_txt(oc_txt) if os.path.isfile(oc_txt) else ('', '')
        parts = []
        if k_locus or k_type:
            parts.append(f"{k_locus}({k_type})" if k_locus and k_type else k_locus or k_type)
        if oc_locus or oc_type:
            parts.append(f"{oc_locus}({oc_type})" if oc_locus and oc_type else oc_locus or oc_type)
        return '/'.join(parts) if parts else 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse Kaptive output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_vibrio_serotype(sample_dir, sample_id):
    k_txt = os.path.join(sample_dir, 'kaptive_vp', f'{sample_id}_vp_k.txt')
    o_txt = os.path.join(sample_dir, 'kaptive_vp', f'{sample_id}_vp_o.txt')
    if not os.path.isfile(k_txt) and not os.path.isfile(o_txt):
        return None
    try:
        k_locus, k_type = _parse_kaptive_txt(k_txt) if os.path.isfile(k_txt) else ('', '')
        o_locus, o_type = _parse_kaptive_txt(o_txt) if os.path.isfile(o_txt) else ('', '')
        parts = []
        if k_locus or k_type:
            parts.append(f"{k_locus}({k_type})" if k_locus and k_type else k_locus or k_type)
        if o_locus or o_type:
            parts.append(f"{o_locus}({o_type})" if o_locus and o_type else o_locus or o_type)
        return '/'.join(parts) if parts else 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse Kaptive VP output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_pseudomonas_serotype(sample_dir, sample_id):
    pasty_tsv = os.path.join(sample_dir, 'pasty', f'{sample_id}.tsv')
    if not os.path.isfile(pasty_tsv):
        return None
    try:
        with open(pasty_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                stype = row.get('type', '').strip()
                return stype if stype else 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse pasty TSV for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


def get_listeria_serotype(sample_dir, sample_id):
    liss_txt = os.path.join(sample_dir, 'lissero', 'lissero_output.txt')
    if not os.path.isfile(liss_txt):
        return None
    try:
        with open(liss_txt) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                serotype = row.get('SEROTYPE', '').strip()
                return serotype if serotype else 'Not detected'
    except Exception as e:
        print(f"Warning: Could not parse LisSero output for {sample_id}: {e}", file=sys.stderr)
        return 'Not detected'


# BMGAP2 data parser

def parse_bmgap2(sample_dir, sample_id, scheme, hinfluenzae_txt=None):
    d = {k: NO_DATA for k in [
        'penA_allele', 'penA_mutations', 'penA_phenotype',
        'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
        'parC_allele', 'parC_phenotype',
        'rpoB_allele', 'rpoB_phenotype',
        'ponA_allele', 'ponA_phenotype',
        'predicted_resistance',
        'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc',
        'FHbp_variant', 'FHbp_subfamily', 'FHbp_peptide',
        'NadA_variant', 'NhbA_peptide', 'vaccine_4CMenB_coverage',
        'folA_allele', 'folA_phenotype',
        'blaTEM1_status', 'blaROB1_status',
    ]}

    # AMR JSON
    amr_jsons = glob.glob(os.path.join(sample_dir, 'bmgap2_amr', f'{sample_id}*amr_data.json'))
    if amr_jsons:
        try:
            with open(amr_jsons[0]) as f:
                amr = json.load(f)
            genes = amr.get('amr_genes', {})

            pena = find_gene(genes, 'penA')
            if pena:
                d['penA_allele']    = pena.get('allele', 'Not detected')
                muts = list(pena.get('known_mutations', {}).keys())
                d['penA_mutations'] = ';'.join(muts) if muts else 'None'
                pen_pheno = (amr.get('antimicrobics', {})
                                .get('Penicillins', {}).get('Penicillin', {}))
                d['penA_phenotype'] = pen_pheno.get('predicted_phenotype', 'Not detected')

            gyra = find_gene(genes, 'gyrA')
            if gyra:
                d['gyrA_allele']    = gyra.get('allele', 'Not detected')
                muts = list(gyra.get('known_mutations', {}).keys())
                d['gyrA_mutations'] = ';'.join(muts) if muts else 'None'
                d['gyrA_phenotype'] = 'Susceptible' if not muts else 'Check'

            parc = find_gene(genes, 'parC')
            if parc:
                d['parC_allele']    = parc.get('allele', 'Not detected')
                muts = list(parc.get('known_mutations', {}).keys())
                d['parC_phenotype'] = 'Susceptible' if not muts else 'Check'

            rpob = find_gene(genes, 'rpoB')
            if rpob:
                d['rpoB_allele']    = rpob.get('allele', 'Not detected')
                muts = list(rpob.get('known_mutations', {}).keys())
                d['rpoB_phenotype'] = 'Susceptible' if not muts else 'Resistant'

            pona = find_gene(genes, 'ponA')
            if pona:
                d['ponA_allele']    = pona.get('allele', 'Not detected')
                muts = list(pona.get('known_mutations', {}).keys())
                d['ponA_phenotype'] = 'Susceptible' if not muts else 'Check'

            if scheme == 'hinfluenzae':
                ftsi = find_gene(genes, 'ftsI')
                if ftsi:
                    d['penA_allele']    = ftsi.get('allele', 'Not detected')
                    muts = list(ftsi.get('known_mutations', {}).keys())
                    d['penA_mutations'] = ';'.join(muts) if muts else 'None'
                    pen_dict = amr.get('antimicrobics', {}).get('Penicillins', {})
                    if pen_dict:
                        d['penA_phenotype'] = next(iter(pen_dict.values())).get(
                            'predicted_phenotype', 'Not detected'
                        )
                    else:
                        d['penA_phenotype'] = 'Susceptible' if not muts else 'Check'

                fola = find_gene(genes, 'folA')
                if fola:
                    d['folA_allele']    = fola.get('allele', 'Not detected')
                    muts = list(fola.get('known_mutations', {}).keys())
                    d['folA_phenotype'] = 'Susceptible' if not muts else 'Resistant'

                blatem = find_gene(genes, 'blaTEM-1')
                if blatem:
                    d['blaTEM1_status'] = blatem.get('status', 'Not detected')

                blarob = find_gene(genes, 'blaROB-1')
                if blarob:
                    d['blaROB1_status'] = blarob.get('status', 'Not detected')

            d['predicted_resistance'] = amr.get('summary', {}).get('predicted_resistance', 'None')

        except Exception as e:
            print(f"Warning: Could not parse BMGAP2 AMR JSON for {sample_id}: {e}", file=sys.stderr)

    # LocusExtractor CSV
    le_dirs = glob.glob(os.path.join(sample_dir, 'bmgap2_locusextractor', f'LE_*_{sample_id}_*'))
    if le_dirs:
        le_csv = os.path.join(le_dirs[0], 'Results_text', f'molecular_data_{sample_id}.csv')
        if os.path.isfile(le_csv):
            try:
                with open(le_csv) as f:
                    reader = csv.DictReader(f)
                    all_rows = list(reader)

                prokka_rows = [r for r in all_rows if 'prokka' in r.get('Filename', '')]
                row = prokka_rows[0] if prokka_rows else (all_rows[0] if all_rows else None)
                if row is not None:
                    if scheme == 'hinfluenzae':
                        d['bmgap2_mlst_st'] = row.get('Hi_MLST_ST', '') or NO_DATA
                        hi_st = d['bmgap2_mlst_st']
                        d['bmgap2_mlst_cc'] = NO_DATA
                        if (hinfluenzae_txt and os.path.isfile(hinfluenzae_txt)
                                and hi_st not in [NO_DATA, 'Not detected', 'New', 'NA', '']):
                            try:
                                with open(hinfluenzae_txt) as ht:
                                    for hrow in ht:
                                        hcols = hrow.strip().split('\t')
                                        if hcols[0] == hi_st:
                                            d['bmgap2_mlst_cc'] = (
                                                hcols[8] if len(hcols) >= 9 else 'Not detected'
                                            )
                                            break
                            except Exception:
                                pass
                    else:
                        d['bmgap2_mlst_st'] = row.get('Nm_MLST_ST', '') or NO_DATA
                        d['bmgap2_mlst_cc'] = row.get('Nm_MLST_cc', '') or NO_DATA

                    d['FHbp_variant']  = normalize_le_value(row.get('FHbp_protein_subvariant_Novartis', ''))
                    d['FHbp_subfamily'] = normalize_le_value(row.get('FHbp_subfamily', ''))
                    d['FHbp_peptide']  = normalize_le_value(row.get('FHbp_protein_subvariant_Oxford', ''))
                    d['NadA_variant']  = normalize_le_value(row.get('NadA_Protein_subvariant_Novartis', ''))
                    d['NhbA_peptide']  = normalize_le_value(row.get('NhbA_Protein_subvariant_Novartis', ''))

                    if scheme == 'hinfluenzae':
                        d['vaccine_4CMenB_coverage'] = 'Not applicable'
                    else:
                        has_fhbp = d['FHbp_variant'] not in [NO_DATA, 'Not detected']
                        has_nada = d['NadA_variant'] not in [NO_DATA, 'Not detected']
                        has_nhba = d['NhbA_peptide'] not in [NO_DATA, 'Not detected']
                        if has_fhbp and (has_nada or has_nhba):
                            d['vaccine_4CMenB_coverage'] = 'Likely'
                        elif has_fhbp:
                            d['vaccine_4CMenB_coverage'] = 'Possible'
                        else:
                            d['vaccine_4CMenB_coverage'] = 'Unlikely'
            except Exception as e:
                print(f"Warning: Could not parse LocusExtractor CSV for {sample_id}: {e}", file=sys.stderr)

    # BMScan JSON
    bmscan_dir = os.path.join(sample_dir, 'bmgap2_bmscan')
    if os.path.isdir(bmscan_dir):
        bmscan_jsons = glob.glob(os.path.join(bmscan_dir, 'species_analysis_*.json'))
        if bmscan_jsons:
            try:
                with open(bmscan_jsons[0]) as f:
                    bmscan = json.load(f)
                for sample_data in bmscan.values():
                    if 'mash_results' in sample_data:
                        d['bmgap2_species'] = sample_data['mash_results'].get('species', '-')
                        break
            except Exception as e:
                print(f"Warning: Could not parse BMScan JSON for {sample_id}: {e}", file=sys.stderr)

    return d


# Report headers

HEADER_STANDARD = [
    'sampleID',
    'mash_species', 'mash_reference', 'mash_distance',
    'kraken_species', 'kraken_percent',
    'blast_16s_tophit', 'blast_16s_pident',
    'skani_species', 'skani_ani', 'skani_align_fraction', 'skani_reference',
    'contamination_flag',
    'mlst_scheme', 'mlst_st', 'serotype',
    'num_clean_reads', 'avg_read_length', 'avg_read_qual', 'est_coverage',
    'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length', 'gc_content',
    'annotated_cds',
]

HEADER_NM = [
    'sampleID',
    'pmga_species', 'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc', 'serotype_notes', 'nm_serogroup', 'predicted_resistance',
    'penA_allele', 'penA_mutations', 'penA_phenotype',
    'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
    'parC_allele', 'parC_phenotype',
    'rpoB_allele', 'rpoB_phenotype',
    'ponA_allele', 'ponA_phenotype',
    'FHbp_variant', 'FHbp_subfamily', 'FHbp_peptide',
    'NadA_variant', 'NhbA_peptide', 'vaccine_4CMenB_coverage',
]

HEADER_HI = [
    'sampleID',
    'pmga_species', 'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc', 'serotype_notes', 'hi_serotype', 'predicted_resistance',
    'ftsI_allele', 'ftsI_mutations', 'ftsI_phenotype',
    'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
    'parC_allele', 'parC_phenotype',
    'rpoB_allele', 'rpoB_phenotype',
    'folA_allele', 'folA_phenotype',
    'blaTEM1_status', 'blaROB1_status',
]


# Main

def main():
    parser = argparse.ArgumentParser(description='Build Sanibel summary reports.')
    parser.add_argument('--outdir',          required=True, help='Pipeline output directory')
    parser.add_argument('--neisseria_txt',   default=None,  help='Neisseria MLST CC lookup table')
    parser.add_argument('--hinfluenzae_txt', default=None,  help='H. influenzae MLST CC lookup table')
    args = parser.parse_args()

    outdir          = args.outdir
    neisseria_txt   = args.neisseria_txt
    hinfluenzae_txt = args.hinfluenzae_txt

    samples = sorted(
        f.replace('_assembly_stats.txt', '')
        for f in glob.glob('*_assembly_stats.txt')
    )

    if not samples:
        print('summary_report.py: no *_assembly_stats.txt files found.', file=sys.stderr)
        sys.exit(1)

    rows_std = []
    rows_nm  = []
    rows_hi  = []

    for sid in samples:
        sample_dir = os.path.join(outdir, sid)

        asm  = parse_assembly_stats(f'{sid}_assembly_stats.txt')
        rm   = parse_read_metrics(f'{sid}_readMetrics.txt')
        cds  = parse_prokka_txt(f'{sid}.txt')
        mlst = parse_mlst(f'{sid}.mlst', neisseria_txt, hinfluenzae_txt)
        kr   = parse_kraken_report(f'{sid}.report')
        pmga = parse_pmga(f'{sid}sta.txt')

        agg_path = f'{sid}_candidate_species.txt'
        agg = parse_aggregate_species(agg_path) if os.path.isfile(agg_path) else \
              {'genus': NO_DATA, 'contamination_flag': NO_DATA}

        _sk_empty = {'ani': NO_DATA, 'confirmed_species': NO_DATA, 'align_fraction': NO_DATA, 'reference': NO_DATA}
        skani_path = f'{sid}_skani.tsv'
        sk = parse_skani(skani_path) if os.path.isfile(skani_path) else _sk_empty

        skani_ID_val     = sk['confirmed_species']
        skani_ANI_val    = sk['ani']
        skani_align_val  = sk['align_fraction']
        skani_ID_ref_val = sk['reference']

        skani_genus  = genus_of(skani_ID_val) if skani_ID_val not in (NO_DATA, 'NO ID', '', None) else None
        mash_genus   = asm['genus']
        kraken_genus = kr['species'].split()[0] if kr['species'] != NO_DATA else None

        blast16s_path   = f'{sid}_16s_blast.tsv'
        anchor_16s      = [skani_genus] if skani_genus else [mash_genus, kraken_genus]
        blast16s_result = parse_blast16s_result(blast16s_path, anchor_genera=anchor_16s) \
                          if os.path.isfile(blast16s_path) \
                          else {'pident': NO_DATA, 'tophit': NO_DATA}

        consensus_genus = mash_genus if (mash_genus and kraken_genus and
                                         mash_genus.lower() == kraken_genus.lower()) else None
        if skani_genus and consensus_genus and skani_genus.lower() != consensus_genus.lower() \
                and frozenset({skani_genus.lower(), consensus_genus.lower()}) in SYNONYMOUS_PAIRS:
            skani_ID_val = f"{asm['genus']}_{asm['species']}"

        if skani_genus and os.path.isfile(blast16s_path):
            contamination_flag = detect_contamination(
                extract_contam_candidates(blast16s_path), skani_genus)
        else:
            contamination_flag = agg['contamination_flag']

        scheme  = mlst['scheme']
        pmga_sp = pmga['species'] or NO_DATA

        common = [
            sid,
            f"{asm['genus']}_{asm['species']}", asm['accession'], asm['mash_distance'],
            kr['species'], kr['percent'],
        ]

        if scheme == 'neisseria':
            bm = parse_bmgap2(sample_dir, sid, scheme, hinfluenzae_txt)
            serotype = pmga['prediction'] or NO_DATA
            std_row = common + [
                blast16s_result['tophit'], blast16s_result['pident'],
                skani_ID_val, skani_ANI_val, skani_align_val, skani_ID_ref_val,
                contamination_flag,
                scheme, mlst['st'], serotype,
                rm['num_reads'], rm['avg_read_len'], rm['avg_qual'], rm['coverage'],
                asm['num_contigs'], asm['longest_contig'], asm['n50'], asm['l50'],
                asm['total_length'], asm['gc_content'], cds,
            ]
            rows_std.append(std_row)
            nm_row = [
                sid,
                pmga_sp,
                bm['bmgap2_species'], bm['bmgap2_mlst_st'], bm['bmgap2_mlst_cc'],
                pmga['serotype_notes'],
                pmga['prediction'] or NO_DATA,
                bm['predicted_resistance'],
                bm['penA_allele'], bm['penA_mutations'], bm['penA_phenotype'],
                bm['gyrA_allele'], bm['gyrA_mutations'], bm['gyrA_phenotype'],
                bm['parC_allele'], bm['parC_phenotype'],
                bm['rpoB_allele'], bm['rpoB_phenotype'],
                bm['ponA_allele'], bm['ponA_phenotype'],
                bm['FHbp_variant'], bm['FHbp_subfamily'], bm['FHbp_peptide'],
                bm['NadA_variant'], bm['NhbA_peptide'], bm['vaccine_4CMenB_coverage'],
            ]
            rows_nm.append(nm_row)

        elif scheme == 'hinfluenzae':
            bm = parse_bmgap2(sample_dir, sid, scheme, hinfluenzae_txt)
            serotype = pmga['prediction'] or NO_DATA
            std_row = common + [
                blast16s_result['tophit'], blast16s_result['pident'],
                skani_ID_val, skani_ANI_val, skani_align_val, skani_ID_ref_val,
                contamination_flag,
                scheme, mlst['st'], serotype,
                rm['num_reads'], rm['avg_read_len'], rm['avg_qual'], rm['coverage'],
                asm['num_contigs'], asm['longest_contig'], asm['n50'], asm['l50'],
                asm['total_length'], asm['gc_content'], cds,
            ]
            rows_std.append(std_row)
            hi_row = [
                sid,
                pmga_sp,
                bm['bmgap2_species'], bm['bmgap2_mlst_st'], bm['bmgap2_mlst_cc'],
                pmga['serotype_notes'],
                pmga['prediction'] or NO_DATA,
                bm['predicted_resistance'],
                bm['penA_allele'], bm['penA_mutations'], bm['penA_phenotype'],
                bm['gyrA_allele'], bm['gyrA_mutations'], bm['gyrA_phenotype'],
                bm['parC_allele'], bm['parC_phenotype'],
                bm['rpoB_allele'], bm['rpoB_phenotype'],
                bm['folA_allele'], bm['folA_phenotype'],
                bm['blaTEM1_status'], bm['blaROB1_status'],
            ]
            rows_hi.append(hi_row)

        else:
            # Species-specific serotype from published output dirs
            serotype = NO_DATA
            for getter in [
                lambda: get_ecoli_serotype(sample_dir, sid),
                lambda: get_klebsiella_serotype(sample_dir, sid),
                lambda: get_legionella_serotype(sample_dir, sid),
                lambda: get_salmonella_serotype(sample_dir, sid),
                lambda: get_gas_serotype(sample_dir, sid),
                lambda: get_shigella_serotype(sample_dir, sid),
                lambda: get_pneumococcal_serotype(sample_dir, sid),
                lambda: get_acinetobacter_serotype(sample_dir, sid),
                lambda: get_vibrio_serotype(sample_dir, sid),
                lambda: get_pseudomonas_serotype(sample_dir, sid),
                lambda: get_listeria_serotype(sample_dir, sid),
            ]:
                result = getter()
                if result is not None:
                    serotype = result
                    break

            row = common + [
                blast16s_result['tophit'], blast16s_result['pident'],
                skani_ID_val, skani_ANI_val, skani_align_val, skani_ID_ref_val,
                contamination_flag,
                scheme, mlst['st'], serotype,
                rm['num_reads'], rm['avg_read_len'], rm['avg_qual'], rm['coverage'],
                asm['num_contigs'], asm['longest_contig'], asm['n50'], asm['l50'],
                asm['total_length'], asm['gc_content'], cds,
            ]
            rows_std.append(row)

    def write_report(path, header, rows):
        with open(path, 'w') as fh:
            fh.write('\t'.join(header) + '\n')
            for row in sorted(rows, key=lambda r: r[0]):
                fh.write('\t'.join(str(v) for v in row) + '\n')
        print(f"summary_report.py: wrote {path} ({len(rows)} sample(s))")

    if rows_std:
        write_report('sum_report.txt',    HEADER_STANDARD, rows_std)
    if rows_nm:
        write_report('nm_sum_report.txt', HEADER_NM,       rows_nm)
    if rows_hi:
        write_report('hi_sum_report.txt', HEADER_HI,       rows_hi)

    if not (rows_std or rows_nm or rows_hi):
        print('summary_report.py: no rows generated.', file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
