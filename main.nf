#!/usr/bin/env nextflow

process echo{
    echo true
    output:
    stdout to out
    
    script:
    """
    samtools --version
    
    gatk --help
    """
}
