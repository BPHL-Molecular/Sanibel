process bmgap2_amr {
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
    import gzip
    import shutil
    import json 

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
        print(f"BMGAP2-AMR: Error reading pyoutputs - {str(e)}", file=sys.stderr)
        sys.exit(0)  # Exit gracefully
    
    if scheme not in ['neisseria', 'hinfluenzae']:
        print(f"BMGAP2-AMR: Skipping {sample_name} - not a meningitis species (scheme={scheme})")
        sys.exit(0)
    
    species_map = {
        'neisseria': 'Nm',
        'hinfluenzae': 'Hi'
    }
    species_code = species_map[scheme]
    
    pmga_dir = f"${mypath}/pmga"
    amr_output_dir = f"${mypath}/bmgap2_amr"
    meningitis_scripts_dir = "${params.bmgap2_db}"
    runast_script = os.path.join(meningitis_scripts_dir, "amr_variants_component", "runAST.py")
    
    if not os.path.isdir(pmga_dir):
        print(f"BMGAP2-AMR: Warning - PMGA directory not found: {pmga_dir}")
        print(f"BMGAP2-AMR: Skipping AMR analysis for {sample_name}")
        sys.exit(0)
     
    json_files = glob.glob(f"{pmga_dir}/*blast-final-results.json*")
    if not json_files:
        print(f"BMGAP2-AMR: Warning - No PMGA JSON file found in {pmga_dir}")
        print(f"BMGAP2-AMR: Skipping AMR analysis for {sample_name}")
        sys.exit(0)
    
    pmga_json_file = json_files[0]
    print(f"BMGAP2-AMR: Found PMGA JSON: {pmga_json_file}")
    
    pmga_json_dir = f"{pmga_dir}/json"
    os.makedirs(pmga_json_dir, exist_ok=True)
    
    original_basename = os.path.basename(pmga_json_file).replace('.gz', '')
    fixed_basename = original_basename.replace('-final-results.json', '_final_results.json')
    json_target = os.path.join(pmga_json_dir, fixed_basename)
    
    if not os.path.exists(json_target):
        if pmga_json_file.endswith('.gz'):
            print(f"BMGAP2-AMR: Decompressing {pmga_json_file} to {fixed_basename}")
            with gzip.open(pmga_json_file, 'rb') as f_in:
                with open(json_target, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        else:
            os.symlink(pmga_json_file, json_target)
    
    with open(json_target, 'r') as f:
        json_data = json.load(f)
    
    needs_transform = False
    
    if isinstance(json_data, dict):
        if 'species' in json_data and 'contigs' in json_data:
            needs_transform = True
            print(f"BMGAP2-AMR: Detected Sanibel format (top-level species/contigs)")
        else:
            first_key = list(json_data.keys())[0]
            if isinstance(json_data[first_key], dict) and 'contigs' in json_data[first_key]:
                needs_transform = False
                print(f"BMGAP2-AMR: Detected BMGAP2 format (sample ID wrapper)")
            else:
                needs_transform = True
                print(f"BMGAP2-AMR: Detected old Sanibel format (flat contigs)")
    
    if needs_transform:
        print(f"BMGAP2-AMR: Transforming to BMGAP2 format")
        
        species_names = {
            'Nm': 'neisseria',
            'Hi': 'hinfluenzae'
        }
        
        if 'species' in json_data and 'contigs' in json_data:
            transformed = {
                sample_name: json_data
            }
            contig_count = len(json_data['contigs'])
        else:
            transformed = {
                sample_name: {
                    'contigs': json_data,
                    'species': species_names.get(species_code, 'unknown')
                }
            }
            contig_count = len(json_data)
        
        with open(json_target, 'w') as f:
            json.dump(transformed, f, indent=2)
        
        print(f"BMGAP2-AMR: Transformation complete - {contig_count} contigs wrapped")
    else:
        print(f"BMGAP2-AMR: JSON already in BMGAP2 format, skipping transformation")
    
    os.makedirs(amr_output_dir, exist_ok=True)
    
    if not os.path.isfile(runast_script):
        print(f"BMGAP2-AMR: Error - runAST.py not found at {runast_script}", file=sys.stderr)
        print(f"BMGAP2-AMR: Check params.bmgap2_db setting: {meningitis_scripts_dir}")
        sys.exit(0)  # Exit gracefully

    print(f"BMGAP2-AMR: Running AMR mutation analysis for {sample_name} (species: {species_code})")
    
    cmd = [
        sys.executable,  # Use same Python as Nextflow (with biopython)
        runast_script,
        "-i", pmga_json_dir,  # Pass json subdirectory
        "-o", amr_output_dir,
        "-s", species_code
    ]
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        
        amr_json = os.path.join(amr_output_dir, f"{sample_name}_amr_data.json")
        if os.path.isfile(amr_json):
            print(f"BMGAP2-AMR: Successfully created {amr_json}")
        else:
            print(f"BMGAP2-AMR: Warning - Expected output file not found: {amr_json}")

    except subprocess.CalledProcessError as e:
        print(f"BMGAP2-AMR: Error running runAST.py:", file=sys.stderr)
        print(e.stderr, file=sys.stderr)
        print(f"BMGAP2-AMR: Continuing pipeline despite AMR analysis failure")
        sys.exit(0)  # Exit gracefully

    except Exception as e:
        print(f"BMGAP2-AMR: Unexpected error - {str(e)}", file=sys.stderr)
        sys.exit(0)  # Exit gracefully
    
    """
}
