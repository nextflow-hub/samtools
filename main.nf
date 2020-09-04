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
        .set { ch_in_samtools }

/*
#==============================================
samtools
#==============================================
*/

process samtoolsFaidx {
    publishDir params.resultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'


    input:
    path refFasta from ch_refFasta

    output:
    file('*.fai') into ch_out_samtools

    script:

    """
    samtools faidx $params.refFasta
    """
}


/*
#==============================================
# extra
#==============================================
*/
