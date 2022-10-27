process pyTask2 {
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
    filepath = "${mypath}/"+items[-1]+"_readMetrics.txt"
    with open(filepath, 'r') as metrics:
        firstline = metrics.readline()
        firstline = firstline.rstrip()
        firstline = firstline.split()
        secondline = metrics.readline()
        secondline = secondline.rstrip()
        secondline = secondline.split()
        rm = dict(zip(firstline, secondline))
        avg_read_len = rm['avgReadLength']
        avg_qual = rm['avgQuality']
        num_read = rm['numReads']
        cov = rm['coverage']
    #print(cov)
    
    #f = open("pyoutputs2.txt", "w")
    f = open("${pyoutputs}", "a")
    f.write(","+str(avg_read_len)+","+str(avg_qual)+","+str(num_read)+","+str(cov))
    f.close
    """
}
