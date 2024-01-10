#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PANGENE } from './workflows/pangene.nf'

workflow {
    PFR_PANGENE()
}

workflow PFR_PANGENE {
    PANGENE()
}
