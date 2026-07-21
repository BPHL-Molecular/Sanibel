#!/usr/bin/env python3
"""
run_bmgap2_bmscan.py - species identification for Neisseria / H. influenzae samples.

Reads the MLST scheme from --mlst to determine whether to run.
Exits cleanly (code 0) if the sample is not a meningitis species.
"""

import argparse
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from bmgap2_helpers import read_scheme, run_tool

TAG = "BMGAP2-BMScan"


def main():
    parser = argparse.ArgumentParser(description='BMGAP2 BMScan species identification.')
    parser.add_argument('--sample',       required=True, help='Sample ID')
    parser.add_argument('--mlst',         required=True, help='MLST result file')
    parser.add_argument('--assembly-dir', required=True, help='Directory holding the assembly FASTAs')
    parser.add_argument('--outdir',       required=True, help='Directory to write BMScan results into')
    parser.add_argument('--db',           required=True, help='BMGAP2 analysis_scripts directory')
    parser.add_argument('--cpus',         required=True, help='Threads for BMScan')
    args = parser.parse_args()

    sample_name = args.sample
    read_scheme(args.mlst, sample_name, TAG)

    assembly_dir = os.path.abspath(args.assembly_dir)
    output_dir   = os.path.abspath(args.outdir)
    bmscan       = os.path.join(args.db, "SpeciesDB", "bin", "identify_species.py")

    if not os.path.isfile(bmscan):
        print(f"{TAG}: Error - BMScan script not found: {bmscan}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)

    print(f"{TAG}: Running species identification for {sample_name}")

    cmd = [
        sys.executable,
        bmscan,
        "-d", assembly_dir,
        "-o", output_dir,
        "-t", args.cpus,
        "-j",
    ]

    run_tool(cmd, TAG, on_fail='fail')

    # BMScan names its JSON by run, so prefix it or samples collide when staged together.
    target = os.path.join(output_dir, f"{sample_name}_species_analysis.json")
    produced = [p for p in glob.glob(os.path.join(output_dir, 'species_analysis_*.json'))
                if p != target]
    if not produced:
        print(f"{TAG}: Error - No species_analysis JSON produced for {sample_name}", file=sys.stderr)
        sys.exit(1)

    shutil.move(produced[0], target)
    print(f"{TAG}: Successfully completed species identification for {sample_name}")


if __name__ == "__main__":
    main()
