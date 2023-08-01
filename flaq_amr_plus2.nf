#!/usr/bin/env nextflow

/*
Note:
Before running the script, please set the parameters in the config file params.yaml
*/

//Step1:input data files
nextflow.enable.dsl=2
def L001R1Lst = []
def sampleNames = []
myDir = file("$params.input")

myDir.eachFileMatch ~/.*_1.fastq.gz/, {L001R1Lst << it.name}
L001R1Lst.sort()
L001R1Lst.each{
   def x = it.minus("_1.fastq.gz")
     //println x
   sampleNames.add(x)
}
//println L001R1Lst
//println sampleNames


//Step2: process the inputed data
A = Channel.fromList(sampleNames)
//A.view()

//include { assemble } from './modules/assemble.nf'
include { fastqc } from './modules/fastqc.nf'
include { trimmomatic } from './modules/trimmomatic.nf'
include { bbtools } from './modules/bbtools.nf'
include { fastqc2 } from './modules/fastqc2.nf'
include { multiqc } from './modules/multiqc.nf'
include { mash } from './modules/mash.nf'
include { unicycler } from './modules/unicycler.nf'
//include { pmga } from './modules/pmga.nf'
include { quast } from './modules/quast.nf'

include { pyTask1 } from './modules/pyTask1.nf'
include { readssum } from './modules/readssum.nf'
include { pyTask2 } from './modules/pyTask2.nf'

include { prokka } from './modules/prokka.nf'
include { amrfinder } from './modules/amrfinder.nf'
include { mlst } from './modules/mlst.nf'
include { kraken } from './modules/kraken.nf'

include { pyTask3 } from './modules/pyTask3.nf'
include { pyTask4 } from './modules/pyTask4.nf'
include { plusAnalyses } from './modules/plusAnalyses.nf'

workflow {
    fastqc(A) | trimmomatic | bbtools| fastqc2 | multiqc | mash | unicycler | quast | pyTask1 | readssum | pyTask2 | prokka | amrfinder | mlst | kraken | pyTask3 | pyTask4 | plusAnalyses | view
}
