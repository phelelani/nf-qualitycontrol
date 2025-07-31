#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================================

// PROCESSES
process run_fastp {
    label 'fastp'
    tag { sample_id }
    publishDir "${qcoutdir}/fastp/${sample_id}", mode: 'copy', overwrite: false
    
    input:
    tuple val(sample_id), path(reads), val(qcoutdir)
    
    output:
    path("*.{html,json}"), emit: fastqc_report
    
    """
    /home/phelelani/applications/fastp/fastp \
        --in1 ${reads.get(0)} --in2 ${reads.get(1)} \
        --json "${sample_id}_fastp.json" --html "${sample_id}_fastp.html" \
        --report_title "${sample_id}_QC_Report" \
        --thread ${task.cpus}
    """
}

// =====================================================================================================
// END
