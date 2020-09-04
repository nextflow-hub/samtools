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


params.index = false
params.faidx = false
params.sort = false

params.bwaMemResultsDir = 'results/bwa/mem'
params.samtoolsFaidxResultsDir = 'results/samtools/faidx'
params.samtoolsSortResultsDir = 'results/samtools/sort'
params.samtoolsIndexResultsDir = 'results/samtools/index'

params.saveMode = 'copy'

params.refFasta = "NC000962_3.fasta"
params.readsFilePattern = "./*_{R1,R2}.fastq.gz"
params.bamFilePattern = ".bam"
params.sortedBamFilePattern = ".sort.bam"


Channel.value("$workflow.launchDir/$params.refFasta")
        .set { ch_refFasta }

Channel.fromFilePairs(params.readsFilePattern)
        .set { ch_in_samtools }

Channel.fromPath("${params.bwaMemResultsDir}/*${params.bamFilePattern}")
        .set { ch_in_samtoolsSort }

Channel.fromPath("${params.samtoolsSortResultsDir}/*${params.sortedBamFilePattern}")
        .set { ch_in_samtoolsIndex }


/*
#==============================================
samtools
#==============================================
*/

process samtoolsFaidx {
    publishDir params.samtoolsFaidxResultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'

    when:
    params.faidx

    input:
    path refFasta from ch_refFasta

    output:
    file('*.fai') into ch_out_samtoolsFaidx

    script:

    """
    samtools faidx $params.refFasta
    """
}



process samtoolsSort {
    publishDir params.samtoolsSortResultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'

    when:
    params.sort

    input:
    file(bamRead) from ch_in_samtoolsSort

    output:
    file("*${params.sortedBamFilePattern}") into ch_out_samtoolsSort

    script:

    genomeName = bamRead.toString().split("\\.")[0]
    """
    samtools sort ${bamRead} >  ${genomeName}.sort.bam
    """
}





process samtoolsIndex {
    publishDir params.samtoolsIndexResultsDir, mode: params.saveMode
    container 'quay.io/biocontainers/samtools:1.10--h2e538c0_3'

    when:
    params.index

    input:
    path refFasta from ch_refFasta
    file(sortedBam) from ch_in_samtoolsIndex

    output:
    file("*.bai") into ch_out_samtoolsIndex

    script:

    """
    samtools index ${sortedBam}
    """
}

/*
#==============================================
# extra
#==============================================
*/

