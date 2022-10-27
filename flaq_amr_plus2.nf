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

include { assemble } from './modules/assemble.nf'
include { pyTask1 } from './modules/pyTask1.nf'
include { readssum } from './modules/readssum.nf'
include { pyTask2 } from './modules/pyTask2.nf'
include { annotate } from './modules/annotate.nf'
include { pyTask3 } from './modules/pyTask3.nf'
include { pyTask4 } from './modules/pyTask4.nf'
include { plusAnalyses } from './modules/plusAnalyses.nf'

workflow {
    assemble(A) | pyTask1 | readssum | pyTask2 |  annotate| pyTask3 | pyTask4 | plusAnalyses | view
}
