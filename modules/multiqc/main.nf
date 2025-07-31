#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================================

// PROCESSES
process run_multiqc {
    label 'multiqc'
    tag { 'multiqc:annotation' }
    publishDir "${qcoutdir}/${tool}", mode: 'copy', overwrite: false
    
    input:
    tuple val(tool), val(qcoutdir), path(data)

    output:
    path("${tool}_combined"), emit: multiqc_report

    """
    multiqc . --force --outdir ${tool}_combined
    """
}

// =====================================================================================================
// END

