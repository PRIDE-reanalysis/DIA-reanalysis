#!/usr/bin/env nextflow

/*
===============================================================================
            PRIDE re-SWATH downstream workflow v1
===============================================================================
 @authors
 Mathias Walzer <walzer@ebi.ac.uk>
-------------------------------------------------------------------------------
Pipeline overview:
 - 0:   MSstats process
 - 1:   Normalisation and exports
 - 2:   `quality` report
-------------------------------------------------------------------------------
Some notes on the technical aspects of the workflow:
 - make sure the container has all required R-libs installed
*/

def helpMessage() {
    log.info"""
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run downstream.nf -c nextflow.config -profile local

    Workflow arguments:

    --annotation The dataset annotation file for MSstats
    --tric The result file from TRIC
    --descriptor The descriptor for the dataset appearing in file names and the report PDF
    --out_dir The destination folder (abs) for the report and exports

    """.stripIndent()
}


/*
 * Config Variables
 */
// show help message if requested
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// defaults where conditionaly needed
params.descriptor = "PXD123"

// prepare input (into value & channels)
annot_channel_1 = Channel.fromPath( "${params.annotation}" )
annot_channel_2 = Channel.fromPath( "${params.annotation}" )
tric_channel = Channel.fromPath( "${params.tric}" )

/*
 * Process definitions
 */
process tricMSSTATS {
    container "${params.aRgh.container}"
    memory { 12.GB * task.attempt } 
    errorStrategy 'retry'
    maxRetries 3

    input:
    file(tric) from tric_channel
    file(annotation) from annot_channel_1

    output:
    file '*.rda' into rda_channel

    script:
    """
    Rscript --vanilla /scripts/DIA_downstream_process.R -a ${annotation} -t ${tric}
    """
}


process rdaMSSTATS {
    container "${params.aRgh.container}"
    memory { 8.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3
    publishDir "${params.out_dir}/" , mode: 'copy', pattern: "*.txt"
    publishDir "${params.out_dir}/" , mode: 'copy', pattern: "*.tsv"
    publishDir "${params.out_dir}/" , mode: 'copy', pattern: "*.pdf"

    input:
    file(rda_file) from rda_channel
    file(annotation) from annot_channel_2

    output:
    file('*.txt') into pg_channel
    file('*.tsv') into tab_channel
    file('*.pdf') into report_channel

    script:
    """
    # N.B.: the directory of wherever this rmd is in will be the working directory by default! e.g. /tmp/this.rmd will operate in /tmp/ and _not_ from where the Rcmd is being called from 
    cp /scripts/DIA_downstream_report.rmd .
    python3 /scripts/rq_pxd.py -a ${params.descriptor}.meta
    Rscript --vanilla -e "rmarkdown::render('DIA_downstream_report.rmd', params = list(rda = '${rda_file}', ann = '${annotation}', idf = '${params.descriptor}.meta', pxdid = '${params.descriptor}') )"
    """
    // no direct knit from /script/DIA_downstream_report.rmd - read-only filesystem
    // also means all output is being routed to /tmp/
    // works but copies rmd (big)
    //sed 's/REPLACEME/${params.descriptor}/' /scripts/diamsstats_qc.rmd > diamsstats_qc_PXD.rmd
    //Rscript --vanilla -e "rmarkdown::render('diamsstats_qc_PXD.rmd', params = list(rda = '\"\$(readlink -f ${rda_file})\"', ann = '\$(readlink -f ${annotation})\' , out = 'proteingroups.txt') )"
}

