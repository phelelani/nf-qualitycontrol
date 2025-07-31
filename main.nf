#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ====================================================================================================
//
//                                         DEFAULT PATAMETERS
//
// ====================================================================================================
workflow            = null
project             = null
outdir              = null
samplesheet_fastqs  = null
// ====================================================================================================


// ====================================================================================================
//
//        SELECT WHICH WORKFLOW TO RUN - USES DIFFERENT INPUT, BUT ALL START FROM SAMPLE SHEET!
//
// ====================================================================================================
workflow                  = params.workflow
project                   = params.project

outdir                    = file(params.outdir, type: 'dir')
samplesheetdir            = file(params.outdir + '/' + params.project + '/_samplesheets', type: 'dir')
qcoutdir                  = file(params.outdir + '/' + params.project + '/_qc_reports', type: 'dir')

outdir.mkdirs()
samplesheetdir.mkdirs()
qcoutdir.mkdirs()
// ====================================================================================================


// ====================================================================================================
//
//        SELECT WHICH WORKFLOW TO RUN - USES DIFFERENT INPUT, BUT ALL START FROM SAMPLE SHEET!
//
// ====================================================================================================
//


// ====================================================================================================
//
//        SELECT WHICH WORKFLOW TO RUN - USES DIFFERENT INPUT, BUT ALL START FROM SAMPLE SHEET!
//
// ====================================================================================================
include { run_fastqc } from '../nf-modules/fastqc/main.nf'
include { run_fastp } from '../nf-modules/fastp/main.nf'
include { run_multiqc } from '../nf-modules/multiqc/main.nf'
// =====================================================================================================


// ====================================================================================================
//
//        SELECT WHICH WORKFLOW TO RUN - USES DIFFERENT INPUT, BUT ALL START FROM SAMPLE SHEET!
//
// ====================================================================================================
// FASTQC WORKFLOW
// ****************************************************************************************************
workflow FASTQC {
    take:
    samples_qc

    main:
    samples_qc
        .map { it -> [ it[0], it[1], qcoutdir ] }
        .set { samples_qc_input }
    run_fastqc(samples_qc_input)

    Channel.of( [ 'fastqc', qcoutdir ])
        .concat (run_fastqc.out.fastqc_report.collect().toSortedList())
        .collect()
        .set { fastqc_multiqc_files }
    run_multiqc(fastqc_multiqc_files)
}

// ****************************************************************************************************
// FASTP WORKFLOW
// ****************************************************************************************************
workflow FASTP {
    take:
    samples_qc

    main:
    samples_qc
        .map { it -> [ it[0], it[1], qcoutdir ] }
        .set { samples_qc_input }
    run_fastp(samples_qc_input)

    Channel.of( [ 'fastp', qcoutdir ])
        .concat (run_fastp.out.fastqc_report.collect().toSortedList())
        .collect()
        .set { fastp_multiqc_files }
    run_multiqc(fastp_multiqc_files)
}

// ====================================================================================================
//
//        SELECT WHICH WORKFLOW TO RUN - USES DIFFERENT INPUT, BUT ALL START FROM SAMPLE SHEET!
//
// ====================================================================================================
workflow {
    switch (workflow) {
            // ****************************************************************************************
            // RUN FASTQC WORKFLOW
            // ****************************************************************************************
        case ['fastqc']:
            // GET SAMPLESHEET            
            samplesheet = Channel.fromPath(params.samplesheet_fastqs, checkIfExists: true) 
            samplesheet
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.Gender}", "${row.FastqR1}", "${row.FastqR2}",
                                "${row.Flowcell}", "${row.Lane}", "${row.BAM}", "${row.gVCF}" ] }
                .set { samples_info }

            samples_info
                .map { [ it[0], [ it[2], it[3] ] ] }
                .set { samples_qc }

            // RUN WORKFLOW
            FASTQC(samples_qc)

            // PRINT SUMMARY
            // workflowSummary()
            break
            // ****************************************************************************************
            // RUN FASTP WORKFLOW
            // ****************************************************************************************
        case ['fastp']:
            // GET SAMPLESHEET            
            samplesheet = Channel.fromPath(params.samplesheet_fastqs, checkIfExists: true) 
            samplesheet
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.Gender}", "${row.FastqR1}", "${row.FastqR2}",
                                "${row.Flowcell}", "${row.Lane}", "${row.BAM}", "${row.gVCF}" ] }
                .set { samples_info }

            samples_info
                .map { [ it[0], [ it[2], it[3] ] ] }
                .set { samples_qc }

            // RUN WORKFLOW
            FASTP(samples_qc)

            // PRINT SUMMARY
            // workflowSummary()
            break
            // ****************************************************************************************
            // RUN EVERYTHIN
            // ****************************************************************************************            
        case ['all']:
            // GET SAMPLESHEET
            samplesheet = Channel.fromPath(params.samplesheet_fastqs, checkIfExists: true)
            samplesheet
                .splitCsv(header: true, sep: '\t')
                .map { row -> [ "${row.SampleID}", "${row.Gender}", "${row.FastqR1}", "${row.FastqR2}",
                                "${row.Flowcell}", "${row.Lane}", "${row.BAM}", "${row.gVCF}" ] }
                .set { samples_info }

            // INPUT FOR FASTQC
            samples_info
                .map { [ it[0], [ it[2], it[3] ] ] }
                .set { samples_qc }

            // INPUT FOR BWAMEM
            samples_info
                .map { [ it[0], it[2], it[3], it[4], it[5] ] }
                .set { samples_bwamem }

            // RUN ALL THE WORKFLOWS
            FASTQC(samples_qc)
            FASTP(samples_qc)       
            // PRINT SUMMARY
            // workflowSummary()
            break
            // ****************************************************************************************
            // ANYTHING ELESE!
            // ****************************************************************************************
        default:
            exit 1, "NO WORKFLOW GIVEN!"
            break
            // ****************************************************************************************
    }
}
// ====================================================================================================
