process bmgap2_locusextractor {
    input:
        val mypath
        path pyoutputs
    output:
        val mypath
        path pyoutputs

    when:
        params.meningitis

    """
    #!/usr/bin/env python3

    import subprocess
    import os
    import sys
    import glob

    items = "${mypath}".strip().split("/")
    sample_name = items[-1]

    pyoutputs_file = "${pyoutputs}"
    scheme = ""

    try:
        with open(pyoutputs_file, 'r') as f:
            content = f.read().strip()
            fields = content.split(',')
            if len(fields) > 16:
                scheme = fields[16]
    except Exception as e:
        print(f"BMGAP2-LocusExtractor: Error reading pyoutputs - {str(e)}", file=sys.stderr)
        sys.exit(0)

    if scheme not in ['neisseria', 'hinfluenzae']:
        print(f"BMGAP2-LocusExtractor: Skipping {sample_name} - not a meningitis species (scheme={scheme})")
        sys.exit(0)

    assembly_dir = f"${mypath}/{sample_name}_assembly"
    output_dir = f"${mypath}/bmgap2_locusextractor"
    meningitis_scripts_dir = "${params.bmgap2_db}"
    locusextractor_dir = os.path.join(meningitis_scripts_dir, "locusextractor")
    locusextractor_script = os.path.join(locusextractor_dir, "LocusExtractor_RHEL8.py")
    settings_dir = os.path.join(locusextractor_dir, "settings")

    if not os.path.isdir(assembly_dir):
        print(f"BMGAP2-LocusExtractor: Error - Assembly directory not found: {assembly_dir}", file=sys.stderr)
        sys.exit(0)

    if not os.path.isfile(locusextractor_script):
        print(f"BMGAP2-LocusExtractor: Error - LocusExtractor script not found: {locusextractor_script}", file=sys.stderr)
        sys.exit(0)

    os.makedirs(output_dir, exist_ok=True)

    print(f"BMGAP2-LocusExtractor: Running vaccine antigen typing for {sample_name}")

    cmd = [
        sys.executable,  # Use same Python as Nextflow (with biopython)
        locusextractor_script,
        "--no_update",  # Use local references, don't download
        "-s", settings_dir,
        "-p", sample_name,
        assembly_dir  # Pass directory, LocusExtractor will scan for FASTA files
    ]

    try:
        result = subprocess.run(
            cmd,
            cwd=output_dir,
            capture_output=True,
            text=True,
            check=False
        )

        if result.returncode != 0:
            print(f"BMGAP2-LocusExtractor: Warning - LocusExtractor returned error code {result.returncode}", file=sys.stderr)
            print(f"STDERR: {result.stderr}", file=sys.stderr)

        print(result.stdout)

        le_dirs = glob.glob(os.path.join(output_dir, f"LE_*_{sample_name}_*"))
        if le_dirs:
            le_dir = le_dirs[0]  # Get most recent (should only be one)
            expected_output = os.path.join(le_dir, "Results_text", f"molecular_data_{sample_name}.csv")
            if os.path.isfile(expected_output):
                print(f"BMGAP2-LocusExtractor: Successfully generated {expected_output}")
            else:
                print(f"BMGAP2-LocusExtractor: Warning - Expected output not found: {expected_output}", file=sys.stderr)
        else:
            print(f"BMGAP2-LocusExtractor: Warning - LocusExtractor output directory not found", file=sys.stderr)

    except Exception as e:
        print(f"BMGAP2-LocusExtractor: Error running LocusExtractor - {str(e)}", file=sys.stderr)
        sys.exit(0)

    """
}
