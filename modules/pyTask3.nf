process pyTask3 {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    """
    #!/usr/bin/env python3

    import subprocess
    items = "${mypath}".strip().split("/")
    filepath1="${mypath}/"+items[-1]+"_assembly/prokka/"+items[-1]+".txt"
    filepath2="${mypath}/"+items[-1]+"_assembly/"+items[-1]+".mlst"
    filepath3="${mypath}/kraken_out/"+items[-1]+".report"
    filepath4="${mypath}/report.txt"
    uppath="/".join(items[:-2])
    filepath5a=uppath+"/neisseria.txt"
    filepath5b=uppath+"/hinfluenzae.txt"
    filepath6="${mypath}/pmga/"+items[-1]+"sta.txt"
    
    cds = ''
    with open(filepath1, 'r') as genes:
        for line in genes:
            line = line.rstrip()
            content = line.split()
            if content[0] == 'CDS:':
                cds = content[1]
    with open(filepath2,'r') as mlst:
       for line in mlst:
           out = line.rstrip().split()
           scheme = out[1]
           st = out[2]
           cc = ""
           if scheme == "neisseria":
              with open(filepath5a, 'r') as mlst_table:
                  for row in mlst_table:
                      cols = row.rstrip().split()
                      if st == cols[0]:
                          cc = cols[8]
                          break

              #serotype prediction of Neisseria 
              subprocess.run('singularity exec -B ' + "${mypath}"+':/data docker://staphb/pmga:latest pmga --blastdir /pmga/blastdbs -o /data/pmga --force --species neisseria /data/'+items[-1]+'_assembly/'+items[-1]+'.fasta', shell=True, check=True) 
                     
           if scheme == "hinfluenzae":
              with open(filepath5b, 'r') as mlst_table:
                  for row in mlst_table:
                      cols = row.rstrip().split()
                      if st == cols[0]:
                          cc = cols[8]
                          break
              #serotype prediction of Hinfluenzae
              subprocess.run('singularity exec -B ' + "${mypath}"+':/data docker://staphb/pmga:latest pmga --blastdir /pmga/blastdbs -o /data/pmga --force --species hinfluenzae /data/'+items[-1]+'_assembly/'+items[-1]+'.fasta', shell=True, check=True) 

    with open(filepath3, 'r') as kreport:
        lines = kreport.readlines()
        for l in lines:
            l_parse = l.lstrip().rstrip().split("\t")
            percent = l_parse[0]
            tax_level = l_parse[3]
            tax = l_parse[5].lstrip()
            if tax_level == 'S':
                break
        
    import os.path
    pgspecies=''
    pgpredict=''
    if os.path.exists(filepath6):
        with open(filepath6, 'r') as pgmalines:
            pglines = pgmalines.readlines()
            pgcells = pglines[1].strip().split("\t")
            pgpredict = pgcells[2]
            pgspecies = pgcells[1]
                 
    #print(cds+"\t"+st+"\t"+tax)
    
    f = open("${pyoutputs}", "a")
    f.write(","+str(cds)+","+str(scheme)+","+str(st)+","+str(cc)+","+str(percent)+","+str(tax)+","+str(pgspecies)+","+str(pgpredict))
    f.close
    
    #f=open("${pyoutputs}", "r")
    #print(f.read())
    
    """
}
