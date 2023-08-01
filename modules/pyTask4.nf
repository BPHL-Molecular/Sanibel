process pyTask4 {
    input:
        val mypath
        path pyoutputs
    output:
        //stdout
        val mypath
        path pyoutputs
        
    $/
    #!/usr/bin/env python3
    items = "${mypath}".strip().split("/")
    #parentpath = "/".join(items[:-1])
    
    with open("${pyoutputs}", "r") as aline:
        for line in aline:
            cells=line.rstrip().split(",")
            results=[items[-1],cells[0]+'_'+cells[1],cells[3],cells[2],cells[20],cells[19],cells[16],cells[17],cells[18],cells[21],cells[22],cells[13],cells[11],cells[12],cells[14],cells[5],cells[6],cells[7],cells[8],cells[9],cells[10],cells[15]]
            print(results)

    report = open("${mypath}"+"/report.txt", 'w')
    header = ['sampleID', 'speciesID_mash', 'nearest_neighbor_mash', 'mash_distance', 'speciesID_kraken', 'kraken_percent', 'mlst_scheme', 'mlst_st', 'mlst_cc', 'pmga_species', 'serotype', 'num_clean_reads', 'avg_readlength', 'avg_read_qual', 'est_coverage', 'num_contigs', 'longest_contig', 'N50', 'L50', 'total_length', 'gc_content', 'annotated_cds'] 
    report.write("\t".join(header))
    report.write('\n')
    report.write('\t'.join(results))
    report.close()
    /$
}
