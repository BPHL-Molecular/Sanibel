process pyTask4 {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    $/
    #!/usr/bin/env python3
    import json
    import os
    import sys
    import glob
    import csv

    items = "${mypath}".strip().split("/")
    sample_id = items[-1]

    with open("${pyoutputs}", "r") as aline:
        for line in aline:
            cells = line.rstrip().split(",")
            results = [
                items[-1],                    # sampleID
                cells[13],                    # num_clean_reads
                cells[11],                    # avg_readlength
                cells[12],                    # avg_read_qual
                cells[14],                    # est_coverage
                cells[5],                     # num_contigs
                cells[6],                     # longest_contig
                cells[7],                     # N50
                cells[8],                     # L50
                cells[9],                     # total_length
                cells[10],                    # gc_content
                cells[15],                    # annotated_cds
                cells[0]+'_'+cells[1],        # speciesID_mash
                cells[3],                     # nearest_neighbor_mash
                cells[2],                     # mash_distance
                cells[20],                    # speciesID_kraken
                cells[19],                    # kraken_percent
                cells[16],                    # mlst_scheme
                cells[17],                    # mlst_st
                cells[18],                    # mlst_cc
                cells[21],                    # pmga_species
                cells[22],                    # serotype
            ]
            scheme = cells[16]

    serotype_notes = '-'
    pmga_sta_file = os.path.join("${mypath}", "pmga", f"{sample_id}sta.txt")
    if os.path.isfile(pmga_sta_file):
        try:
            with open(pmga_sta_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    data_line = lines[1].rstrip().split('\t')
                    if len(data_line) >= 5:
                        serotype_notes = data_line[4]
        except Exception as e:
            print(f"Warning: Could not parse PMGA sta file for {sample_id}: {e}", file=sys.stderr)

    results.append(serotype_notes)

    # Sentinel values:
    # 'No data'        = module did not run or output file missing
    # 'Not detected'   = module ran, gene/variant not found
    # 'Not applicable' = field does not apply to this species
    bmgap2_data = {
        'penA_allele': 'No data',
        'penA_mutations': 'No data',
        'penA_phenotype': 'No data',
        'gyrA_allele': 'No data',
        'gyrA_mutations': 'No data',
        'gyrA_phenotype': 'No data',
        'parC_allele': 'No data',
        'parC_phenotype': 'No data',
        'rpoB_allele': 'No data',
        'rpoB_phenotype': 'No data',
        'ponA_allele': 'No data',
        'ponA_phenotype': 'No data',
        'predicted_resistance': 'No data',
        'bmgap2_species': 'No data',
        'bmgap2_mlst_st': 'No data',
        'bmgap2_mlst_cc': 'No data',
        'FHbp_variant': 'No data',
        'FHbp_subfamily': 'No data',
        'FHbp_peptide': 'No data',
        'NadA_variant': 'No data',
        'NhbA_peptide': 'No data',
        'vaccine_4CMenB_coverage': 'No data',
        'folA_allele': 'No data',      # Hi-specific
        'folA_phenotype': 'No data',   # Hi-specific
        'blaTEM1_status': 'No data',   # Hi-specific
        'blaROB1_status': 'No data'    # Hi-specific
    }

    if "${params.meningitis}" == "true" and scheme in ['neisseria', 'hinfluenzae']:

        amr_base = os.path.join("${mypath}", "bmgap2_amr")
        amr_json_files = glob.glob(os.path.join(amr_base, f"{sample_id}*amr_data.json"))

        if amr_json_files:
            try:
                with open(amr_json_files[0], 'r') as f:
                    amr_data = json.load(f)

                def find_gene(gene_dict, gene_name):
                    for key, value in gene_dict.items():
                        if value.get('Gene_name') == gene_name:
                            return value
                    for key, value in gene_dict.items():
                        if key == gene_name or f"({gene_name})" in key or key.startswith(gene_name + " "):
                            return value
                    return None

                def normalize_le_value(val):
                    if not val:
                        return 'Not detected'
                    if val in ['Not found', 'not found']:
                        return 'Not detected'
                    if val.startswith('Allele not identified') or val.startswith('Peptide not found'):
                        return 'Not detected'
                    return val

                amr_genes = amr_data.get('amr_genes', {})

                pena = find_gene(amr_genes, 'penA')
                if pena:
                    bmgap2_data['penA_allele'] = pena.get('allele', 'Not detected')
                    known_muts = list(pena.get('known_mutations', {}).keys())
                    bmgap2_data['penA_mutations'] = ';'.join(known_muts) if known_muts else 'None'
                    if 'Penicillin' in amr_data.get('antimicrobics', {}).get('Penicillins', {}):
                        bmgap2_data['penA_phenotype'] = amr_data['antimicrobics']['Penicillins']['Penicillin'].get('predicted_phenotype', 'Not detected')

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

                # ftsI reuses penA_* slots (Hi only — Hi has no penA)
                if scheme == 'hinfluenzae':
                    ftsi = find_gene(amr_genes, 'ftsI')
                    if ftsi:
                        bmgap2_data['penA_allele'] = ftsi.get('allele', 'Not detected')
                        known_muts = list(ftsi.get('known_mutations', {}).keys())
                        bmgap2_data['penA_mutations'] = ';'.join(known_muts) if known_muts else 'None'
                        if 'Penicillins' in amr_data.get('antimicrobics', {}):
                            for drug, drug_data in amr_data['antimicrobics']['Penicillins'].items():
                                bmgap2_data['penA_phenotype'] = drug_data.get('predicted_phenotype', 'Not detected')
                                break
                        else:
                            bmgap2_data['penA_phenotype'] = 'Susceptible' if not known_muts else 'Check'

                summary = amr_data.get('summary', {})
                bmgap2_data['predicted_resistance'] = summary.get('predicted_resistance', 'None')

                # Hi-specific AMR genes
                if scheme == 'hinfluenzae':
                    fola = find_gene(amr_genes, 'folA')
                    if fola:
                        bmgap2_data['folA_allele'] = fola.get('allele', 'Not detected')
                        known_muts = list(fola.get('known_mutations', {}).keys())
                        bmgap2_data['folA_phenotype'] = 'Susceptible' if not known_muts else 'Resistant'

                    blatem = find_gene(amr_genes, 'blaTEM-1')
                    if blatem:
                        bmgap2_data['blaTEM1_status'] = blatem.get('status', 'Not detected')

                    blarob = find_gene(amr_genes, 'blaROB-1')
                    if blarob:
                        bmgap2_data['blaROB1_status'] = blarob.get('status', 'Not detected')

            except Exception as e:
                print(f"Warning: Could not parse AMR JSON for {sample_id}: {e}", file=sys.stderr)

        le_base = os.path.join("${mypath}", "bmgap2_locusextractor")
        le_dirs = glob.glob(os.path.join(le_base, f"LE_*_{sample_id}_*"))
        le_csv_path = os.path.join(le_dirs[0], "Results_text", f"molecular_data_{sample_id}.csv") if le_dirs else None

        if le_csv_path and os.path.isfile(le_csv_path):
            try:
                with open(le_csv_path, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        if 'Assembly-1' in row.get('Unique_ID', ''):
                            if scheme == 'hinfluenzae':
                                bmgap2_data['bmgap2_mlst_st'] = row.get('Hi_MLST_ST', '') or 'No data'
                                hi_st = bmgap2_data['bmgap2_mlst_st']
                                uppath = "/".join(items[:-2])
                                hi_txt = uppath + "/hinfluenzae.txt"
                                bmgap2_data['bmgap2_mlst_cc'] = 'No data'
                                if os.path.isfile(hi_txt) and hi_st not in ['No data', 'Not detected', 'New', 'NA', '']:
                                    try:
                                        with open(hi_txt, 'r') as ht:
                                            for hrow in ht:
                                                hcols = hrow.rstrip().split('\t')
                                                if hcols[0] == hi_st:
                                                    bmgap2_data['bmgap2_mlst_cc'] = hcols[8] if len(hcols) >= 9 else 'Not detected'
                                                    break
                                    except Exception:
                                        pass
                            else:
                                bmgap2_data['bmgap2_mlst_st'] = row.get('Nm_MLST_ST', '') or 'No data'
                                bmgap2_data['bmgap2_mlst_cc'] = row.get('Nm_MLST_cc', '') or 'No data'

                            bmgap2_data['FHbp_variant'] = normalize_le_value(row.get('FHbp_protein_subvariant_Novartis', ''))
                            bmgap2_data['FHbp_subfamily'] = normalize_le_value(row.get('FHbp_subfamily', ''))
                            bmgap2_data['FHbp_peptide'] = normalize_le_value(row.get('FHbp_protein_subvariant_Oxford', ''))
                            bmgap2_data['NadA_variant'] = normalize_le_value(row.get('NadA_Protein_subvariant_Novartis', ''))
                            bmgap2_data['NhbA_peptide'] = normalize_le_value(row.get('NhbA_Protein_subvariant_Novartis', ''))

                            if scheme == 'hinfluenzae':
                                bmgap2_data['vaccine_4CMenB_coverage'] = 'Not applicable'
                            else:
                                has_fhbp = bmgap2_data['FHbp_variant'] not in ['No data', 'Not detected']
                                has_nada = bmgap2_data['NadA_variant'] not in ['No data', 'Not detected']
                                has_nhba = bmgap2_data['NhbA_peptide'] not in ['No data', 'Not detected']
                                if has_fhbp and (has_nada or has_nhba):
                                    bmgap2_data['vaccine_4CMenB_coverage'] = 'Likely'
                                elif has_fhbp:
                                    bmgap2_data['vaccine_4CMenB_coverage'] = 'Possible'
                                else:
                                    bmgap2_data['vaccine_4CMenB_coverage'] = 'Unlikely'
                            break

            except Exception as e:
                print(f"Warning: Could not parse LocusExtractor CSV for {sample_id}: {e}", file=sys.stderr)

        bmscan_base = os.path.join("${mypath}", "bmgap2_bmscan")
        if os.path.isdir(bmscan_base):
            bmscan_json_files = glob.glob(os.path.join(bmscan_base, "species_analysis_*.json"))
            if bmscan_json_files:
                try:
                    with open(bmscan_json_files[0], 'r') as f:
                        bmscan_data = json.load(f)
                    for sample_key, sample_data in bmscan_data.items():
                        if 'mash_results' in sample_data:
                            bmgap2_data['bmgap2_species'] = sample_data['mash_results'].get('species', '-')
                            break
                except Exception as e:
                    print(f"Warning: Could not parse BMScan JSON for {sample_id}: {e}", file=sys.stderr)

        results.extend([
            bmgap2_data['bmgap2_species'],
            bmgap2_data['bmgap2_mlst_st'],
            bmgap2_data['bmgap2_mlst_cc'],
            bmgap2_data['predicted_resistance'],
            bmgap2_data['penA_allele'],
            bmgap2_data['penA_mutations'],
            bmgap2_data['penA_phenotype'],
            bmgap2_data['gyrA_allele'],
            bmgap2_data['gyrA_mutations'],
            bmgap2_data['gyrA_phenotype'],
            bmgap2_data['parC_allele'],
            bmgap2_data['parC_phenotype'],
            bmgap2_data['rpoB_allele'],
            bmgap2_data['rpoB_phenotype'],
            bmgap2_data['ponA_allele'],
            bmgap2_data['ponA_phenotype'],
            bmgap2_data['FHbp_variant'],
            bmgap2_data['FHbp_subfamily'],
            bmgap2_data['FHbp_peptide'],
            bmgap2_data['NadA_variant'],
            bmgap2_data['NhbA_peptide'],
            bmgap2_data['vaccine_4CMenB_coverage'],
            bmgap2_data['folA_allele'],
            bmgap2_data['folA_phenotype'],
            bmgap2_data['blaTEM1_status'],
            bmgap2_data['blaROB1_status']
        ])

    # results layout (indices):
    #  0        sampleID
    #  1-11     QC metrics
    #  12-14    speciesID_mash, nearest_neighbor_mash, mash_distance
    #  15-16    speciesID_kraken, kraken_percent
    #  17-18    mlst_scheme, mlst_st
    #  19       mlst_cc
    #  20       pmga_species
    #  21       serotype  (PMGA)
    #  22       serotype_notes
    #  23+      bmgap2 fields (only when meningitis=true and Nm/Hi)

    # Nm and Hi (meningitis=true) go to species-specific reports; all others to report.txt
    is_meningitis_species = "${params.meningitis}" == "true" and scheme in ['neisseria', 'hinfluenzae']

    if not is_meningitis_species:
        gen_header = [
            'sampleID',
            'num_clean_reads', 'avg_readlength', 'avg_read_qual',
            'est_coverage', 'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length',
            'gc_content', 'annotated_cds',
            'speciesID_mash', 'nearest_neighbor_mash', 'mash_distance',
            'speciesID_kraken', 'kraken_percent', 'mlst_scheme', 'mlst_st',
            'serotype'
        ]
        serotype_value = 'No data'

        # E. coli: O:H serotype from SerotypeFinder
        ecoli_json = os.path.join("${mypath}", "escherichia", "data.json")
        if os.path.isfile(ecoli_json):
            try:
                with open(ecoli_json, 'r') as f:
                    ecoli_data = json.load(f)
                sf_results = ecoli_data.get('serotypefinder', {}).get('results', {})
                # Deduplicate: multiple genes (e.g. wzx + wzy) may confirm the same type
                o_types = list(dict.fromkeys(
                    hit.get('serotype', '') for hit in sf_results.get('O_type', {}).values()
                    if hit.get('serotype')
                ))
                h_types = list(dict.fromkeys(
                    hit.get('serotype', '') for hit in sf_results.get('H_type', {}).values()
                    if hit.get('serotype')
                ))
                o_str = '/'.join(sorted(set(o_types))) if o_types else 'NT'
                h_str = '/'.join(sorted(set(h_types))) if h_types else 'NT'
                serotype_value = f"{o_str}:{h_str}"
            except Exception as e:
                print(f"Warning: Could not parse SerotypeFinder JSON for {sample_id}: {e}", file=sys.stderr)
                serotype_value = 'Not detected'

        # Klebsiella: ST from Kleborate
        kleb_tsv = os.path.join("${mypath}", "klebsiella", "kleborate-test-out.tsv")
        if os.path.isfile(kleb_tsv):
            try:
                with open(kleb_tsv, 'r') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        st = row.get('ST', '').strip()
                        serotype_value = st if st and st != '-' else 'Not detected'
                        break
            except Exception as e:
                print(f"Warning: Could not parse Kleborate TSV for {sample_id}: {e}", file=sys.stderr)
                serotype_value = 'Not detected'

        # Legionella: SBT + all 7 loci from Legsta
        legsta_txt = os.path.join("${mypath}", "legsta", "legsta_output.txt")
        if os.path.isfile(legsta_txt):
            try:
                with open(legsta_txt, 'r') as f:
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
                        serotype_value = ','.join(loci)
                        break
            except Exception as e:
                print(f"Warning: Could not parse Legsta output for {sample_id}: {e}", file=sys.stderr)
                serotype_value = 'Not detected'

        # Salmonella: antigenic profile + serotype name from SeqSero2
        salm_dirs = glob.glob(os.path.join("${mypath}", "salmonella", "SeqSero_result_*"))
        if salm_dirs:
            salm_tsv = os.path.join(salm_dirs[0], "SeqSero_result.tsv")
            if os.path.isfile(salm_tsv):
                try:
                    with open(salm_tsv, 'r') as f:
                        reader = csv.DictReader(f, delimiter='\t')
                        for row in reader:
                            profile = row.get('Predicted antigenic profile', '').strip()
                            serotype_name = row.get('Predicted serotype', '').strip()
                            if profile:
                                serotype_value = f"{profile}({serotype_name})" if serotype_name else profile
                            else:
                                serotype_value = 'Not detected'
                            break
                except Exception as e:
                    print(f"Warning: Could not parse SeqSero2 TSV for {sample_id}: {e}", file=sys.stderr)
                    serotype_value = 'Not detected'

        # Streptococcus pyogenes/dysgalactiae: emm type from emm-typing-tool
        gas_txt = os.path.join("${mypath}", "groupAstrep", "groupAstrep_result.txt")
        if os.path.isfile(gas_txt):
            try:
                validated = None
                nonvalidated = None
                with open(gas_txt, 'r') as f:
                    for line in f:
                        fields = line.strip().split()
                        if len(fields) >= 3:
                            emm_raw = fields[2].replace('.sds', '')
                            if fields[1] == 'EMM_validated' and validated is None:
                                validated = emm_raw
                            elif fields[1] == 'EMM_nonValidated' and nonvalidated is None:
                                nonvalidated = emm_raw
                serotype_value = validated or nonvalidated or 'Not detected'
            except Exception as e:
                print(f"Warning: Could not parse emm-typing-tool output for {sample_id}: {e}", file=sys.stderr)
                serotype_value = 'Not detected'

        # Shigella: prediction from Shigatyper
        shigella_txt = os.path.join("${mypath}", "shigella", "shigella_output.txt")
        if os.path.isfile(shigella_txt):
            try:
                with open(shigella_txt, 'r') as f:
                    lines = f.readlines()
                found_header = False
                for line in lines:
                    if found_header:
                        fields = line.strip().split('\t')
                        serotype_value = fields[1].strip() if len(fields) >= 2 else 'Not detected'
                        break
                    if line.startswith('sample\t'):
                        found_header = True
                if not found_header:
                    serotype_value = 'Not detected'
            except Exception as e:
                print(f"Warning: Could not parse Shigatyper output for {sample_id}: {e}", file=sys.stderr)
                serotype_value = 'Not detected'

        # results[0:19] = sampleID through mlst_st
        gen_data = results[0:19] + [serotype_value]
        with open("${mypath}/report.txt", 'w') as report:
            report.write('\t'.join(gen_header) + '\n')
            report.write('\t'.join(gen_data) + '\n')

    # Species-specific meningitis reports
    if is_meningitis_species:
        sp_base_header = [
            'sampleID',
            'num_clean_reads', 'avg_readlength', 'avg_read_qual',
            'est_coverage', 'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length',
            'gc_content', 'annotated_cds',
            'speciesID_mash', 'nearest_neighbor_mash', 'mash_distance',
            'speciesID_kraken', 'kraken_percent', 'mlst_scheme', 'mlst_st', 'mlst_cc',
            'pmga_species'
        ]
        sp_base_data = results[0:21]  # sampleID through pmga_species (indices 0-20)

        if scheme == 'neisseria':
            sp_header = sp_base_header + [
                'nm_serogroup', 'serotype_notes',
                'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc', 'predicted_resistance',
                'penA_allele', 'penA_mutations', 'penA_phenotype',
                'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
                'parC_allele', 'parC_phenotype',
                'rpoB_allele', 'rpoB_phenotype',
                'ponA_allele', 'ponA_phenotype',
                'FHbp_variant', 'FHbp_subfamily', 'FHbp_peptide',
                'NadA_variant', 'NhbA_peptide', 'vaccine_4CMenB_coverage'
            ]
            sp_data = sp_base_data + [
                results[21],   # nm_serogroup
                results[22],   # serotype_notes
                bmgap2_data['bmgap2_species'],
                bmgap2_data['bmgap2_mlst_st'],
                bmgap2_data['bmgap2_mlst_cc'],
                bmgap2_data['predicted_resistance'],
                bmgap2_data['penA_allele'],
                bmgap2_data['penA_mutations'],
                bmgap2_data['penA_phenotype'],
                bmgap2_data['gyrA_allele'],
                bmgap2_data['gyrA_mutations'],
                bmgap2_data['gyrA_phenotype'],
                bmgap2_data['parC_allele'],
                bmgap2_data['parC_phenotype'],
                bmgap2_data['rpoB_allele'],
                bmgap2_data['rpoB_phenotype'],
                bmgap2_data['ponA_allele'],
                bmgap2_data['ponA_phenotype'],
                bmgap2_data['FHbp_variant'],
                bmgap2_data['FHbp_subfamily'],
                bmgap2_data['FHbp_peptide'],
                bmgap2_data['NadA_variant'],
                bmgap2_data['NhbA_peptide'],
                bmgap2_data['vaccine_4CMenB_coverage']
            ]
            with open("${mypath}/report_nm.txt", 'w') as sp_report:
                sp_report.write('\t'.join(sp_header) + '\n')
                sp_report.write('\t'.join(sp_data) + '\n')

        elif scheme == 'hinfluenzae':
            sp_header = sp_base_header + [
                'hi_serotype', 'serotype_notes',
                'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc', 'predicted_resistance',
                'ftsI_allele', 'ftsI_mutations', 'ftsI_phenotype',
                'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
                'parC_allele', 'parC_phenotype',
                'rpoB_allele', 'rpoB_phenotype',
                'folA_allele', 'folA_phenotype',
                'blaTEM1_status', 'blaROB1_status'
            ]
            sp_data = sp_base_data + [
                results[21],   # hi_serotype
                results[22],   # serotype_notes
                bmgap2_data['bmgap2_species'],
                bmgap2_data['bmgap2_mlst_st'],
                bmgap2_data['bmgap2_mlst_cc'],
                bmgap2_data['predicted_resistance'],
                bmgap2_data['penA_allele'],   # ftsI data stored in penA_* slots
                bmgap2_data['penA_mutations'],
                bmgap2_data['penA_phenotype'],
                bmgap2_data['gyrA_allele'],
                bmgap2_data['gyrA_mutations'],
                bmgap2_data['gyrA_phenotype'],
                bmgap2_data['parC_allele'],
                bmgap2_data['parC_phenotype'],
                bmgap2_data['rpoB_allele'],
                bmgap2_data['rpoB_phenotype'],
                bmgap2_data['folA_allele'],
                bmgap2_data['folA_phenotype'],
                bmgap2_data['blaTEM1_status'],
                bmgap2_data['blaROB1_status']
            ]
            with open("${mypath}/report_hi.txt", 'w') as sp_report:
                sp_report.write('\t'.join(sp_header) + '\n')
                sp_report.write('\t'.join(sp_data) + '\n')
    /$
}
