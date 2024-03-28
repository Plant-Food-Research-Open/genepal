#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { validateParameters    } from 'plugin/nf-validation'

validateParameters()

include { PANGENE               } from './workflows/pangene.nf'

workflow {
    PFR_PANGENE()
}

workflow PFR_PANGENE {
    PANGENE()
}
