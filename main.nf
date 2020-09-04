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

params.refFasta = "NC000962_3.fasta"

Channel.value("$workflow.launchDir/$params.refFasta")
        .set { ch_refFasta }

Channel.fromFilePairs(params.filePattern)
        .into { ch_in_samtools }

/*
#==============================================
samtools
#==============================================
*/

process samtoolsIndex {
    publishDir params.resultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/samtools:0.1.19--hfb9b9cc_8'


    input:
    path refFasta from ch_refFasta

    output:
    tuple file('*.fai') into ch_out_samtools

    script:

    """
    samtools index $params.refFasta
    """
}


/*
#==============================================
# extra
#==============================================
*/
