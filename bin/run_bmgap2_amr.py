#!/usr/bin/env python3
"""
run_bmgap2_amr.py - AMR mutation analysis for Neisseria / H. influenzae samples.

Reads the MLST scheme from --mlst to determine whether to run.
Exits cleanly (code 0) if the sample is not a meningitis species.
"""

import argparse
import glob
import gzip
import json
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from bmgap2_helpers import read_scheme, run_tool

TAG = "BMGAP2-AMR"


def main():
    parser = argparse.ArgumentParser(description='BMGAP2 AMR mutation analysis.')
    parser.add_argument('--sample',   required=True, help='Sample ID')
    parser.add_argument('--mlst',     required=True, help='MLST result file')
    parser.add_argument('--pmga-dir', required=True, help='Directory holding the PMGA BLAST JSON')
    parser.add_argument('--outdir',   required=True, help='Directory to write AMR results into')
    parser.add_argument('--db',       required=True, help='BMGAP2 analysis_scripts directory')
    args = parser.parse_args()

    sample_name = args.sample
    scheme      = read_scheme(args.mlst, sample_name, TAG)

    species_code = {'neisseria': 'Nm', 'hinfluenzae': 'Hi'}[scheme]
    runast       = os.path.join(args.db, "amr_variants_component", "runAST.py")

    json_files = glob.glob(os.path.join(args.pmga_dir, '*blast-final-results.json*'))
    if not json_files:
        print(f"{TAG}: Error - No PMGA JSON file found in {args.pmga_dir}", file=sys.stderr)
        sys.exit(1)

    pmga_json_file = json_files[0]
    print(f"{TAG}: Found PMGA JSON: {pmga_json_file}")

    json_stage_dir = os.path.abspath('pmga_json')
    os.makedirs(json_stage_dir, exist_ok=True)

    original_basename = os.path.basename(pmga_json_file).replace('.gz', '')
    fixed_basename    = original_basename.replace('-final-results.json', '_final_results.json')
    json_target       = os.path.join(json_stage_dir, fixed_basename)

    if pmga_json_file.endswith('.gz'):
        print(f"{TAG}: Decompressing {pmga_json_file} to {fixed_basename}")
        with gzip.open(pmga_json_file, 'rb') as f_in:
            with open(json_target, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copyfile(pmga_json_file, json_target)

    with open(json_target) as f:
        json_data = json.load(f)

    needs_transform = False
    if isinstance(json_data, dict):
        if 'species' in json_data and 'contigs' in json_data:
            needs_transform = True
            print(f"{TAG}: Detected Sanibel format (top-level species/contigs)")
        else:
            first_key = list(json_data.keys())[0]
            if isinstance(json_data[first_key], dict) and 'contigs' in json_data[first_key]:
                print(f"{TAG}: Detected BMGAP2 format (sample ID wrapper)")
            else:
                needs_transform = True
                print(f"{TAG}: Detected old Sanibel format (flat contigs)")

    if needs_transform:
        print(f"{TAG}: Transforming to BMGAP2 format")

        if 'species' in json_data and 'contigs' in json_data:
            transformed  = {sample_name: json_data}
            contig_count = len(json_data['contigs'])
        else:
            transformed  = {
                sample_name: {
                    'contigs': json_data,
                    'species': scheme
                }
            }
            contig_count = len(json_data)

        with open(json_target, 'w') as f:
            json.dump(transformed, f, indent=2)

        print(f"{TAG}: Transformation complete - {contig_count} contigs wrapped")
    else:
        print(f"{TAG}: JSON already in BMGAP2 format, skipping transformation")

    amr_out_dir = os.path.abspath(args.outdir)
    os.makedirs(amr_out_dir, exist_ok=True)

    if not os.path.isfile(runast):
        print(f"{TAG}: Error - runAST.py not found at {runast}", file=sys.stderr)
        print(f"{TAG}: Check bmgap2_db setting: {args.db}", file=sys.stderr)
        sys.exit(1)

    print(f"{TAG}: Running AMR mutation analysis for {sample_name} (species: {species_code})")

    cmd = [
        sys.executable,
        runast,
        "-i", json_stage_dir,
        "-o", amr_out_dir,
        "-s", species_code,
    ]

    run_tool(cmd, TAG, on_fail='continue')

    amr_json = os.path.join(amr_out_dir, f"{sample_name}_amr_data.json")
    if os.path.isfile(amr_json):
        print(f"{TAG}: Successfully created {amr_json}")
    else:
        print(f"{TAG}: Warning - Expected output file not found: {amr_json}")


if __name__ == "__main__":
    main()
