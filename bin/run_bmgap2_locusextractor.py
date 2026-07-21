#!/usr/bin/env python3
"""
run_bmgap2_locusextractor.py - vaccine antigen typing for Neisseria / H. influenzae samples.

Reads the MLST scheme from --mlst to determine whether to run.
Exits cleanly (code 0) if the sample is not a meningitis species.
"""

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from bmgap2_helpers import read_scheme, run_tool

TAG = "BMGAP2-LocusExtractor"


def main():
    parser = argparse.ArgumentParser(description='BMGAP2 LocusExtractor vaccine antigen typing.')
    parser.add_argument('--sample',        required=True, help='Sample ID')
    parser.add_argument('--mlst',          required=True, help='MLST result file')
    parser.add_argument('--assembly-dir',  required=True, help='Directory holding the assembly FASTAs')
    parser.add_argument('--outdir',        required=True, help='Directory to write LocusExtractor results into')
    parser.add_argument('--db',            required=True, help='BMGAP2 analysis_scripts directory')
    args = parser.parse_args()

    sample_name = args.sample
    read_scheme(args.mlst, sample_name, TAG)

    assembly_dir = os.path.abspath(args.assembly_dir)
    output_dir   = os.path.abspath(args.outdir)
    le_dir       = os.path.join(args.db, "locusextractor")
    le_script    = os.path.join(le_dir, "LocusExtractor_RHEL8.py")
    settings_dir = os.path.join(le_dir, "settings")

    if not os.path.isfile(le_script):
        print(f"{TAG}: Error - LocusExtractor script not found: {le_script}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    print(f"{TAG}: Running vaccine antigen typing for {sample_name}")

    cmd = [
        sys.executable,
        le_script,
        "--no_update",
        "-s", settings_dir,
        "-p", sample_name,
        assembly_dir,
    ]

    run_tool(cmd, TAG, cwd=output_dir, on_fail='warn')

    le_dirs = glob.glob(os.path.join(output_dir, f"LE_*_{sample_name}_*"))
    if le_dirs:
        expected = os.path.join(
            le_dirs[0], "Results_text", f"molecular_data_{sample_name}.csv"
        )
        if os.path.isfile(expected):
            print(f"{TAG}: Successfully generated {expected}")
        else:
            print(f"{TAG}: Warning - Expected output not found: {expected}", file=sys.stderr)
    else:
        print(f"{TAG}: Warning - LocusExtractor output directory not found", file=sys.stderr)


if __name__ == "__main__":
    main()
