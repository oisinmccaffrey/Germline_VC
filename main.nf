#!/usr/bin/env nextflow

process echo{

    echo true
    output:
    stdout to out
    
    script:
    """
    echo "test Repo"
    """
}
