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
    items = "${mypath}".strip().split("/")
    filepath1="${mypath}/"+items[-1]+"_assembly/prokka/"+items[-1]+".txt"
    filepath2="${mypath}/"+items[-1]+"_assembly/"+items[-1]+".mlst"
    filepath3="${mypath}/kraken_out/"+items[-1]+".report"
    filepath4="${mypath}/report.txt"
    
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
    with open(filepath3, 'r') as kreport:
        lines = kreport.readlines()
        for l in lines:
            l_parse = l.lstrip().rstrip().split("\t")
            percent = l_parse[0]
            tax_level = l_parse[3]
            tax = l_parse[5].lstrip()
            if tax_level == 'S':
                break
    #print(cds+"\t"+st+"\t"+tax)
    
    f = open("${pyoutputs}", "a")
    f.write(","+str(cds)+","+str(scheme)+","+str(st)+","+str(percent)+","+str(tax))
    f.close
    
    #f=open("${pyoutputs}", "r")
    #print(f.read())
    
    """
}
