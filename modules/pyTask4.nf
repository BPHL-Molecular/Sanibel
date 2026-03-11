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
            scheme = cells[16]  # MLST scheme for conditional logic

    serotype_notes = '-'
    pmga_sta_file = os.path.join("${mypath}", "pmga", f"{sample_id}sta.txt")
    if os.path.isfile(pmga_sta_file):
        try:
            with open(pmga_sta_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:  # Skip header, read data line
                    data_line = lines[1].rstrip().split('\t')
                    if len(data_line) >= 5:  # Ensure notes column exists
                        serotype_notes = data_line[4]  # notes is 5th column (index 4)
        except Exception as e:
            print(f"Warning: Could not parse PMGA sta file for {sample_id}: {e}", file=sys.stderr)

    results.append(serotype_notes)

    bmgap2_data = {
        'penA_allele': '-',
        'penA_mutations': '-',
        'penA_phenotype': '-',
        'gyrA_allele': '-',
        'gyrA_mutations': '-',
        'gyrA_phenotype': '-',
        'parC_allele': '-',
        'parC_phenotype': '-',
        'rpoB_allele': '-',
        'rpoB_phenotype': '-',
        'ponA_allele': '-',
        'ponA_phenotype': '-',
        'predicted_resistance': '-',
        'bmgap2_species': '-',
        'bmgap2_mlst_st': '-',
        'bmgap2_mlst_cc': '-',
        'FHbp_variant': '-',
        'FHbp_subfamily': '-',
        'FHbp_peptide': '-',
        'NadA_variant': '-',
        'NhbA_peptide': '-',
        'vaccine_4CMenB_coverage': '-'
    }

    if "${params.meningitis}" == "true" and scheme in ['neisseria', 'hinfluenzae']:

        amr_base = os.path.join("${mypath}", "bmgap2_amr")
        amr_json_files = glob.glob(os.path.join(amr_base, f"{sample_id}*amr_data.json"))

        if amr_json_files:
            amr_json_path = amr_json_files[0]
            try:
                with open(amr_json_path, 'r') as f:
                    amr_data = json.load(f)

                def find_gene(gene_dict, gene_name):
                    for key, value in gene_dict.items():
                        if value.get('Gene_name') == gene_name:
                            return value
                    for key, value in gene_dict.items():
                        if key == gene_name or f"({gene_name})" in key or key.startswith(gene_name + " "):
                            return value
                    return None
                
                amr_genes = amr_data.get('amr_genes', {})

                pena = find_gene(amr_genes, 'penA')
                if pena:
                    bmgap2_data['penA_allele'] = pena.get('allele', '-')
                    known_muts = list(pena.get('known_mutations', {}).keys())
                    bmgap2_data['penA_mutations'] = ';'.join(known_muts) if known_muts else 'None'
                    if 'Penicillin' in amr_data.get('antimicrobics', {}).get('Penicillins', {}):
                        bmgap2_data['penA_phenotype'] = amr_data['antimicrobics']['Penicillins']['Penicillin'].get('predicted_phenotype', '-')

                gyra = find_gene(amr_genes, 'gyrA')
                if gyra:
                    bmgap2_data['gyrA_allele'] = gyra.get('allele', '-')
                    known_muts = list(gyra.get('known_mutations', {}).keys())
                    bmgap2_data['gyrA_mutations'] = ';'.join(known_muts) if known_muts else 'None'
                    bmgap2_data['gyrA_phenotype'] = 'Susceptible' if not known_muts else 'Check'
                
                parc = find_gene(amr_genes, 'parC')
                if parc:
                    bmgap2_data['parC_allele'] = parc.get('allele', '-')
                    known_muts = list(parc.get('known_mutations', {}).keys())
                    bmgap2_data['parC_phenotype'] = 'Susceptible' if not known_muts else 'Check'

                rpob = find_gene(amr_genes, 'rpoB')
                if rpob:
                    bmgap2_data['rpoB_allele'] = rpob.get('allele', '-')
                    known_muts = list(rpob.get('known_mutations', {}).keys())
                    bmgap2_data['rpoB_phenotype'] = 'Susceptible' if not known_muts else 'Resistant'

                pona = find_gene(amr_genes, 'ponA')
                if pona:
                    bmgap2_data['ponA_allele'] = pona.get('allele', '-')
                    known_muts = list(pona.get('known_mutations', {}).keys())
                    bmgap2_data['ponA_phenotype'] = 'Susceptible' if not known_muts else 'Check'

                summary = amr_data.get('summary', {})
                bmgap2_data['predicted_resistance'] = summary.get('predicted_resistance', 'None')
                    
            except Exception as e:
                print(f"Warning: Could not parse AMR JSON for {sample_id}: {e}", file=sys.stderr)

        le_base = os.path.join("${mypath}", "bmgap2_locusextractor")
        le_dirs = glob.glob(os.path.join(le_base, f"LE_*_{sample_id}_*"))
        if le_dirs:
            le_dir = le_dirs[0]  # Get most recent (should only be one)
            le_csv_path = os.path.join(le_dir, "Results_text", f"molecular_data_{sample_id}.csv")
        else:
            le_csv_path = None

        if le_csv_path and os.path.isfile(le_csv_path):
            try:
                import csv
                with open(le_csv_path, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        if 'Assembly-1' in row.get('Unique_ID', ''):
                            bmgap2_data['bmgap2_mlst_st'] = row.get('Nm_MLST_ST', '-')
                            bmgap2_data['bmgap2_mlst_cc'] = row.get('Nm_MLST_cc', '-')

                            bmgap2_data['FHbp_variant'] = row.get('FHbp_protein_subvariant_Novartis', '-')
                            bmgap2_data['FHbp_subfamily'] = row.get('FHbp_subfamily', '-')
                            bmgap2_data['FHbp_peptide'] = row.get('FHbp_protein_subvariant_Oxford', '-')
                            bmgap2_data['NadA_variant'] = row.get('NadA_Protein_subvariant_Novartis', '-')
                            bmgap2_data['NhbA_peptide'] = row.get('NhbA_Protein_subvariant_Novartis', '-')
                            

                            has_fhbp = bmgap2_data['FHbp_variant'] not in ['-', 'Not found', 'Allele not identified']
                            has_nada = bmgap2_data['NadA_variant'] not in ['-', 'Not found']
                            has_nhba = bmgap2_data['NhbA_peptide'] not in ['-', 'Not found']
                            
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
                bmscan_json_path = bmscan_json_files[0]
                try:
                    with open(bmscan_json_path, 'r') as f:
                        bmscan_data = json.load(f)

                    for sample_key, sample_data in bmscan_data.items():
                        if 'mash_results' in sample_data:
                            species = sample_data['mash_results'].get('species', '-')
                            bmgap2_data['bmgap2_species'] = species
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
            bmgap2_data['vaccine_4CMenB_coverage']
        ])
    
    print(results)

    report = open("${mypath}"+"/report.txt", 'w')

    header = [
        'sampleID',
        'num_clean_reads', 'avg_readlength', 'avg_read_qual', 
        'est_coverage', 'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length', 
        'gc_content', 'annotated_cds',
        'speciesID_mash', 'nearest_neighbor_mash', 'mash_distance',
        'speciesID_kraken', 'kraken_percent', 'mlst_scheme', 'mlst_st', 'mlst_cc', 
        'pmga_species', 'serotype', 'serotype_notes'
    ]

    if "${params.meningitis}" == "true":
        header.extend([
            'bmgap2_species', 'bmgap2_mlst_st', 'bmgap2_mlst_cc', 'predicted_resistance',
            'penA_allele', 'penA_mutations', 'penA_phenotype',
            'gyrA_allele', 'gyrA_mutations', 'gyrA_phenotype',
            'parC_allele', 'parC_phenotype',
            'rpoB_allele', 'rpoB_phenotype',
            'ponA_allele', 'ponA_phenotype',
            'FHbp_variant', 'FHbp_subfamily', 'FHbp_peptide',
            'NadA_variant', 'NhbA_peptide', 'vaccine_4CMenB_coverage'
        ])
    report.write("\t".join(header))
    report.write('\n')
    report.write('\t'.join(results))
    report.write('\n')
    report.close()
    /$
}
