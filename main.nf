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
include { SEQKIT } from './modules/nf-core/seqkit'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: Filter Genome Assembly by Minimum Contig Length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Include the seqkit module
include { SEQKIT } from './modules/nf-core/seqkit'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS: Filter Genome Assembly by Minimum Contig Length
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SEQKIT_GET_LENGTH {

    input:
    path input_file

    output:
    path 'filtered_output_file.txt'

    script:
    """
    # Filter contigs based on length and output filtered FASTA
    seqkit seq --min-len ${params.min_contig_length} ${genome_fasta} > filtered_${meta.id}.fasta

    # Generate a list of filtered contigs
    seqkit fx2tab --length --name filtered_${meta.id}.fasta | awk '{print \$1}' > ${meta.id}_contig_list.txt
    """
}

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
    ch_tar_assm_str
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
    // Filter genome assembly by minimum contig length
    //
    SEQKIT_GET_LENGTH(ch_target_assembly)

    //
    // Run GENEPAL main workflow using filtered FASTA
    //
    GENEPAL(
        SEQKIT_GET_LENGTH.out.filtered_fasta.map { meta, fasta, contig_list -> [ meta, fasta ] }, // Filtered genome FASTA
        ch_tar_assm_str,
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
    // SUBWORKFLOW: Run initialization tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.rna_evidence,
        params.liftoff_annotations,
        params.orthofinder_annotations
    )

    //
    // Filter genome assembly by minimum contig length
    //
    SEQKIT_GET_LENGTH(PIPELINE_INITIALISATION.out.target_assembly)

    //
    // Run main workflow using filtered FASTA
    //
    PLANTFOODRESEARCHOPEN_GENEPAL(
        SEQKIT_GET_LENGTH.out.filtered_fasta,
        PIPELINE_INITIALISATION.out.tar_assm_str,
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
    PIPELINE_COMPLETION(
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
