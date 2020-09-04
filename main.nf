#!/usr/bin/env nextflow

/*
#==============================================
code documentation
#==============================================
*/


/*
#==============================================
params
#==============================================
*/

params.resultsDir = 'results/samtools'
params.saveMode = 'copy'
params.filePattern = "./*_{R1,R2}.fastq.gz"


ch_refFILE = Channel.value("$baseDir/refFILE")



Channel.fromFilePairs(params.filePattern)
        .into { ch_in_samtools }

/*
#==============================================
samtools
#==============================================
*/

process samtools {
    publishDir params.resultsDir, mode: params.saveMode
    container 'samtools'


    input:
    set genomeFileName, file(genomeReads) from ch_in_samtools

    output:
    path samtools into ch_out_samtools


    script:
    genomeName = genomeFileName.toString().split("\\_")[0]

    """
    CLI samtools
    """
}


/*
#==============================================
# extra
#==============================================
*/
