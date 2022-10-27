process pyTask1 {
    input:
        val mypath
 
    output:
        //stdout
        val "${mypath}", emit: outputpath2
        path "pyoutputs.txt", emit: pyoutputs
 
    """
    #!/usr/bin/env python3
    import sys
    import os
    import subprocess
    import argparse
    import datetime
    import fnmatch
    import re
    import pandas as pd
    
    items = "${mypath}".strip().split("/")
    #print(items[-1])
    filepath1 = "${mypath}"+"/mash_output/"+items[-1]+"_distances_top10.tab"
    with open(filepath1, 'r') as mash:
        top_hit = mash.readline()
        top_hit = str(top_hit)
        cells = top_hit.split()
        acell = cells[0]
        subcells = re.split('-.-|_subsp',acell)
        subsubcells = subcells[1].split('_')
        genus = subsubcells[0]
        species = subsubcells[1]
        distance = cells[2]
        accession = top_hit.split("-")[5]
    #print(accession)
    
    filepath2 = "${mypath}"+"/"+items[-1]+"_assembly/quast_results/report.tsv"
    df = pd.read_table(filepath2, sep="\t")
    assem = list(df.columns)[1]
    #print(assem)
    contigs = df[assem][12].astype(int)
    long_contig = df[assem][13].astype(int)
    n50 = df[assem][16].astype(int)
    l50 = df[assem][18].astype(int)
    genome = df[assem][14].astype(int)
    gc = df[assem][15].astype(int)
    #print(genome)
    

    f = open("pyoutputs.txt", "w")
    f.write(str(genus)+","+str(species)+","+str(distance)+","+str(accession)+","+str(assem)+","+str(contigs)+","+str(long_contig)+","+str(n50)+","+str(l50)+","+str(genome)+","+str(gc))
    f.close
    
    """
}
