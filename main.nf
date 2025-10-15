#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plant-food-research-open/genepal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/plant-food-research-open/genepal
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENEPAL                   } from './workflows/genepal'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_genepal_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_genepal_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow PLANTFOODRESEARCHOPEN_GENEPAL {

    take:
    ch_target_assembly
    ch_is_masked
    ch_te_library
    ch_braker_annotation
    ch_braker_ex_asm_str
    ch_benchmark_gff
    ch_rna_sra
    ch_rna_fq
    ch_rna_bam_by_assembly
    ch_sortmerna_fastas
    ch_ext_prot_fastas
    ch_liftoff_fasta
    ch_liftoff_gff
    ch_tsebra_config
    ch_orthofinder_pep

    main:
    //
    // WORKFLOW: Run pipeline
    //
    GENEPAL(
        ch_target_assembly,
        ch_is_masked,
        ch_te_library,
        ch_braker_annotation,
        ch_braker_ex_asm_str,
        ch_benchmark_gff,
        ch_rna_sra,
        ch_rna_fq,
        ch_rna_bam_by_assembly,
        ch_sortmerna_fastas,
        ch_ext_prot_fastas,
        ch_liftoff_fasta,
        ch_liftoff_gff,
        ch_tsebra_config,
        ch_orthofinder_pep
    )
    emit:
    multiqc_report = GENEPAL.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.rna_evidence,
        params.liftoff_annotations,
        params.orthofinder_annotations,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    PLANTFOODRESEARCHOPEN_GENEPAL(
        PIPELINE_INITIALISATION.out.target_assembly,
        PIPELINE_INITIALISATION.out.is_masked,
        PIPELINE_INITIALISATION.out.te_library,
        PIPELINE_INITIALISATION.out.braker_annotation,
        PIPELINE_INITIALISATION.out.braker_ex_asm_str,
        PIPELINE_INITIALISATION.out.benchmark_gff,
        PIPELINE_INITIALISATION.out.rna_sra,
        PIPELINE_INITIALISATION.out.rna_fq,
        PIPELINE_INITIALISATION.out.rna_bam_by_assembly,
        PIPELINE_INITIALISATION.out.sortmerna_fastas,
        PIPELINE_INITIALISATION.out.ext_prot_fastas,
        PIPELINE_INITIALISATION.out.liftoff_fasta,
        PIPELINE_INITIALISATION.out.liftoff_gff,
        PIPELINE_INITIALISATION.out.tsebra_config,
        PIPELINE_INITIALISATION.out.orthofinder_pep
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        PLANTFOODRESEARCHOPEN_GENEPAL.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
