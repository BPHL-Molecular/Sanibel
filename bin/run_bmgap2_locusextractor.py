#!/usr/bin/env python3
"""
run_bmgap2_locusextractor.py — vaccine antigen typing for Neisseria / H. influenzae samples.

Usage:
    run_bmgap2_locusextractor.py <mypath> <pyoutputs_file> <bmgap2_db>

Reads MLST scheme from pyoutputs_file to determine whether to run.
Exits cleanly (code 0) if the sample is not a meningitis species.
"""

import glob
import os
import subprocess
import sys


def main():
    if len(sys.argv) != 4:
        print(
            f"Usage: {sys.argv[0]} <mypath> <pyoutputs_file> <bmgap2_db>",
            file=sys.stderr
        )
        sys.exit(1)

    mypath         = sys.argv[1]
    pyoutputs_file = sys.argv[2]
    bmgap2_db      = sys.argv[3]

    items       = mypath.strip().split("/")
    sample_name = items[-1]

    scheme = ""
    try:
        with open(pyoutputs_file) as f:
            content = f.read().strip()
            fields  = content.split(',')
            if len(fields) > 16:
                scheme = fields[16]
    except Exception as e:
        print(f"BMGAP2-LocusExtractor: Error reading pyoutputs - {e}", file=sys.stderr)
        sys.exit(0)

    if scheme not in ['neisseria', 'hinfluenzae']:
        print(
            f"BMGAP2-LocusExtractor: Skipping {sample_name} "
            f"- not a meningitis species (scheme={scheme})"
        )
        sys.exit(0)

    assembly_dir  = f"{mypath}/{sample_name}_assembly"
    output_dir    = f"{mypath}/bmgap2_locusextractor"
    le_dir        = os.path.join(bmgap2_db, "locusextractor")
    le_script     = os.path.join(le_dir, "LocusExtractor_RHEL8.py")
    settings_dir  = os.path.join(le_dir, "settings")

    if not os.path.isdir(assembly_dir):
        print(
            f"BMGAP2-LocusExtractor: Error - Assembly directory not found: {assembly_dir}",
            file=sys.stderr
        )
        sys.exit(0)

    if not os.path.isfile(le_script):
        print(
            f"BMGAP2-LocusExtractor: Error - LocusExtractor script not found: {le_script}",
            file=sys.stderr
        )
        sys.exit(0)

    os.makedirs(output_dir, exist_ok=True)

    print(f"BMGAP2-LocusExtractor: Running vaccine antigen typing for {sample_name}")

    cmd = [
        sys.executable,
        le_script,
        "--no_update",
        "-s", settings_dir,
        "-p", sample_name,
        assembly_dir,
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=output_dir,
            capture_output=True,
            text=True,
            check=False,
        )

        if result.returncode != 0:
            print(
                f"BMGAP2-LocusExtractor: Warning - LocusExtractor returned error code "
                f"{result.returncode}",
                file=sys.stderr
            )
            print(f"STDERR: {result.stderr}", file=sys.stderr)

        print(result.stdout)

        le_dirs = glob.glob(os.path.join(output_dir, f"LE_*_{sample_name}_*"))
        if le_dirs:
            expected = os.path.join(
                le_dirs[0], "Results_text", f"molecular_data_{sample_name}.csv"
            )
            if os.path.isfile(expected):
                print(f"BMGAP2-LocusExtractor: Successfully generated {expected}")
            else:
                print(
                    f"BMGAP2-LocusExtractor: Warning - Expected output not found: {expected}",
                    file=sys.stderr
                )
        else:
            print(
                "BMGAP2-LocusExtractor: Warning - LocusExtractor output directory not found",
                file=sys.stderr
            )

    except Exception as e:
        print(f"BMGAP2-LocusExtractor: Error running LocusExtractor - {e}", file=sys.stderr)
        sys.exit(0)


if __name__ == "__main__":
    main()
