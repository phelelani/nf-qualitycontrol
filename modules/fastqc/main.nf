#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================================

// PROCESSES
process run_fastqc {
    label 'fastqc'
    tag { sample_id }
    publishDir "${qcoutdir}/fastqc/${sample_id}", mode: 'copy', overwrite: false
    
    input:
    tuple val(sample_id), path(reads), val(qcoutdir)
    
    output:
    path("*.{html,zip}"), emit: fastqc_report
    
    """
    fastqc ${reads.findAll().join(' ') } --threads ${task.cpus} --noextract
    """
}

// =====================================================================================================
// END
