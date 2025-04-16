/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                      } from 'plugin/nf-schema'

include { PREPARE_ASSEMBLY                      } from '../subworkflows/local/prepare_assembly'
include { PREPROCESS_RNASEQ                     } from '../subworkflows/local/preprocess_rnaseq'
include { ALIGN_RNASEQ                          } from '../subworkflows/local/align_rnaseq'
include { PREPARE_EXT_PROTS                     } from '../subworkflows/local/prepare_ext_prots'
include { FASTA_BRAKER3                         } from '../subworkflows/local/fasta_braker3'
include { FASTA_LIFTOFF                         } from '../subworkflows/local/fasta_liftoff/main'
include { PURGE_BRAKER_MODELS                   } from '../subworkflows/local/purge_braker_models'
include { GFF_MERGE_CLEANUP                     } from '../subworkflows/local/gff_merge_cleanup'
include { GFF_EGGNOGMAPPER                      } from '../subworkflows/local/gff_eggnogmapper'
include { PURGE_NOHIT_MODELS                    } from '../subworkflows/local/purge_nohit_models'
include { GFF_STORE                             } from '../subworkflows/local/gff_store'
include { FASTA_ORTHOFINDER                     } from '../subworkflows/local/fasta_orthofinder'
include { FASTA_GXF_BUSCO_PLOT                  } from '../subworkflows/nf-core/fasta_gxf_busco_plot/main'

include { GXF_FASTA_AGAT_SPADDINTRONS_SPEXTRACTSEQUENCES                            } from '../subworkflows/gallvp/gxf_fasta_agat_spaddintrons_spextractsequences/main'
include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS as DOWNLOAD_RNASEQ_FROM_SRA  } from '../subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools/main'

include { CAT_CAT as SAVE_MARKED_GFF3           } from '../modules/nf-core/cat/cat/main'
include { GFFCOMPARE as BENCHMARK               } from '../modules/nf-core/gffcompare/main'
include { FILE_GUNZIP as BENCHMARK_GFF3_GUNZIP  } from '../subworkflows/local/file_gunzip'
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { GENEPALREPORT                         } from '../modules/local/genepalreport/main.nf'

include { methodsDescriptionText                } from '../subworkflows/local/utils_nfcore_genepal_pipeline'
include { softwareVersionsToYAML                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { paramsSummaryMultiqc                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENEPAL {

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
    // Channels
    ch_versions                 = Channel.empty()
    ch_multiqc_files            = Channel.empty()

    // SUBWORKFLOW: PREPARE_ASSEMBLY
    PREPARE_ASSEMBLY(
        ch_target_assembly,
        ch_te_library,
        params.repeat_annotator,
        params.repeatmasker_save_outputs,
        ch_braker_ex_asm_str,
        ch_is_masked
    )

    ch_valid_target_assembly    = PREPARE_ASSEMBLY.out.target_assemby
    ch_masked_target_assembly   = PREPARE_ASSEMBLY.out.masked_target_assembly
    ch_target_assemby_index     = PREPARE_ASSEMBLY.out.target_assemby_index
    ch_versions                 = ch_versions.mix(PREPARE_ASSEMBLY.out.versions)

    // SUBWORKFLOW: DOWNLOAD_RNASEQ_FROM_SRA
    DOWNLOAD_RNASEQ_FROM_SRA (ch_rna_sra, [] )

    ch_rna_all_fq               = ch_rna_fq
                                | mix(
                                    DOWNLOAD_RNASEQ_FROM_SRA.out.reads
                                    | map { meta, fq -> [ meta, [ fq ] ] }
                                )
    ch_versions                 = ch_versions.mix(DOWNLOAD_RNASEQ_FROM_SRA.out.versions)

    // SUBWORKFLOW: PREPROCESS_RNASEQ
    PREPROCESS_RNASEQ(
        ch_rna_all_fq,
        ch_tar_assm_str,
        ch_braker_ex_asm_str,
        params.fastqc_skip,
        params.fastp_skip,
        params.save_trimmed,
        params.min_trimmed_reads,
        params.remove_ribo_rna,
        ch_sortmerna_fastas
    )

    ch_trim_reads               = PREPROCESS_RNASEQ.out.trim_reads
    ch_reads_target             = PREPROCESS_RNASEQ.out.reads_target

    ch_multiqc_files            = ch_multiqc_files
                                | mix(PREPROCESS_RNASEQ.out.fastqc_raw_zip)
                                | mix(PREPROCESS_RNASEQ.out.trim_json)
                                | mix(PREPROCESS_RNASEQ.out.fastqc_trim_zip)

    ch_versions                 = ch_versions.mix(PREPROCESS_RNASEQ.out.versions)

    // SUBWORKFLOW: ALIGN_RNASEQ
    ALIGN_RNASEQ(
        ch_reads_target,
        ch_trim_reads,
        ch_rna_bam_by_assembly,
        ch_target_assemby_index,
    )

    ch_rnaseq_bam               = ALIGN_RNASEQ.out.bam

    ch_multiqc_files            = ch_multiqc_files
                                | mix(ALIGN_RNASEQ.out.star_log_final)

    ch_versions                 = ch_versions.mix(ALIGN_RNASEQ.out.versions)

    // MODULE: PREPARE_EXT_PROTS
    PREPARE_EXT_PROTS(
        ch_ext_prot_fastas
    )

    ch_ext_prots_fasta          = PREPARE_EXT_PROTS.out.ext_prots_fasta
    ch_versions                 = ch_versions.mix(PREPARE_EXT_PROTS.out.versions)

    // SUBWORKFLOW: FASTA_BRAKER3
    FASTA_BRAKER3(
        ch_masked_target_assembly,
        ch_braker_ex_asm_str,
        ch_rnaseq_bam,
        ch_ext_prots_fasta,
        ch_braker_annotation
    )

    ch_braker_gff3              = FASTA_BRAKER3.out.braker_gff3
    ch_braker_hints             = FASTA_BRAKER3.out.braker_hints
    ch_versions                 = ch_versions.mix(FASTA_BRAKER3.out.versions)

    // SUBWORKFLOW: FASTA_LIFTOFF
    FASTA_LIFTOFF(
        ch_valid_target_assembly,
        ch_liftoff_fasta,
        ch_liftoff_gff,
        params.filter_liftoff_by_hints,
        ch_braker_hints,
        ch_tsebra_config,
        params.allow_isoforms
    )

    ch_liftoff_gff3             = FASTA_LIFTOFF.out.gff3
    ch_versions                 = ch_versions.mix(FASTA_LIFTOFF.out.versions)

    // SUBWORKFLOW: PURGE_BRAKER_MODELS
    PURGE_BRAKER_MODELS(
        ch_braker_gff3,
        ch_braker_hints,
        ch_liftoff_gff3,
        ch_tsebra_config,
        params.allow_isoforms
    )

    ch_braker_purged_gff        = PURGE_BRAKER_MODELS.out.braker_purged_gff
    ch_versions                 = ch_versions.mix(PURGE_BRAKER_MODELS.out.versions)

    // SUBWORKFLOW: GFF_MERGE_CLEANUP
    GFF_MERGE_CLEANUP(
        ch_braker_purged_gff,
        ch_liftoff_gff3,
        params.filter_genes_by_aa_length,
        params.append_genome_prefix_to_feature_ids
    )

    ch_merged_gff               = GFF_MERGE_CLEANUP.out.gff
    ch_versions                 = ch_versions.mix(GFF_MERGE_CLEANUP.out.versions)

    // SUBWORKFLOW: GFF_EGGNOGMAPPER
    GFF_EGGNOGMAPPER(
        ch_merged_gff,
        ch_valid_target_assembly,
        params.eggnogmapper_db_dir,
    )

    ch_eggnogmapper_hits        = GFF_EGGNOGMAPPER.out.eggnogmapper_hits
    ch_eggnogmapper_annotations = GFF_EGGNOGMAPPER.out.eggnogmapper_annotations
    ch_versions                 = ch_versions.mix(GFF_EGGNOGMAPPER.out.versions)

    // SUBWORKFLOW: PURGE_NOHIT_MODELS
    PURGE_NOHIT_MODELS(
        ch_merged_gff,
        ch_eggnogmapper_hits,
        params.eggnogmapper_purge_nohits && params.eggnogmapper_db_dir
    )

    ch_purged_gff               = PURGE_NOHIT_MODELS.out.purged_gff
    ch_versions                 = ch_versions.mix(PURGE_NOHIT_MODELS.out.versions)

    // SUBWORKFLOW: GFF_STORE
    GFF_STORE(
        ch_purged_gff,
        ch_eggnogmapper_annotations,
        ch_valid_target_assembly,
        params.eggnogmapper_db_dir
    )

    ch_final_gff                = GFF_STORE.out.final_gff
    ch_final_proteins           = GFF_STORE.out.final_proteins
    ch_versions                 = ch_versions.mix(GFF_STORE.out.versions)

    // SUBWORKFLOW: FASTA_ORTHOFINDER
    FASTA_ORTHOFINDER(
        ch_final_proteins,
        ch_orthofinder_pep
    )

    ch_orthofinder_statistics   = FASTA_ORTHOFINDER.out.orthofinder_statistics
    ch_orthofinder_hogs         = FASTA_ORTHOFINDER.out.orthofinder_hogs
    ch_versions                 = ch_versions.mix(FASTA_ORTHOFINDER.out.versions)

    // SUBWORKFLOW: FASTA_GXF_BUSCO_PLOT
    ch_busco_fasta              = params.busco_skip
                                ? Channel.empty()
                                : ch_valid_target_assembly

    ch_busco_gff                = params.busco_skip
                                ? Channel.empty()
                                : ch_final_gff

    FASTA_GXF_BUSCO_PLOT(
        ch_busco_fasta,
        ch_busco_gff,
        'genome',
        params.busco_lineage_datasets?.tokenize(' '),
        [], // val_busco_lineages_path
        [], // val_busco_config
        true
    )

    ch_busco_fasta_summary      = FASTA_GXF_BUSCO_PLOT.out.assembly_short_summaries_txt
    ch_busco_gff_summary        = FASTA_GXF_BUSCO_PLOT.out.annotation_short_summaries_txt

    ch_busco_fasta_plot_summary = FASTA_GXF_BUSCO_PLOT.out.assembly_plot_summary_txt
    ch_busco_gff_plot_summary   = FASTA_GXF_BUSCO_PLOT.out.annotation_plot_summary_txt

    ch_multiqc_files            = ch_multiqc_files
                                | mix(
                                    ch_busco_fasta_plot_summary.map { file -> [ [], file ] }
                                )
                                | mix(
                                    ch_busco_gff_plot_summary.map { file -> [ [], file ] }
                                )

    ch_versions                 = ch_versions.mix(FASTA_GXF_BUSCO_PLOT.out.versions)

    // SUBWORKFLOW: GXF_FASTA_AGAT_SPADDINTRONS_SPEXTRACTSEQUENCES
    GXF_FASTA_AGAT_SPADDINTRONS_SPEXTRACTSEQUENCES(
        ch_final_gff,
        ch_valid_target_assembly
    )

    ch_splicing_marked_gff3     = GXF_FASTA_AGAT_SPADDINTRONS_SPEXTRACTSEQUENCES.out.marked_gff3
    ch_versions                 = ch_versions.mix(GXF_FASTA_AGAT_SPADDINTRONS_SPEXTRACTSEQUENCES.out.versions)

    // MODULE: CAT_CAT as SAVE_MARKED_GFF3
    SAVE_MARKED_GFF3 ( ch_splicing_marked_gff3 )

    ch_saved_marked_gff3        = SAVE_MARKED_GFF3.out.file_out
    ch_versions                 = ch_versions.mix(SAVE_MARKED_GFF3.out.versions.first())


    // SUBWORKFLOW: FILE_GUNZIP as BENCHMARK_GFF3_GUNZIP
    BENCHMARK_GFF3_GUNZIP ( ch_benchmark_gff )
    ch_benchmark_gunzip_gff     = BENCHMARK_GFF3_GUNZIP.out.gunzip
    ch_versions                 = ch_versions.mix(BENCHMARK_GFF3_GUNZIP.out.versions)

    // MODULE: GFFCOMPARE as BENCHMARK
    ch_benchmark_inputs         = ch_final_gff
                                | join ( ch_valid_target_assembly )
                                | join ( ch_benchmark_gunzip_gff )

    BENCHMARK (
        ch_benchmark_inputs.map { meta, gff, fasta, ref_gff -> [ meta, gff ] },
        ch_benchmark_inputs.map { meta, gff, fasta, ref_gff -> [ meta, fasta, [] ] },
        ch_benchmark_inputs.map { meta, gff, fasta, ref_gff -> [ meta, ref_gff ] }
    )

    ch_benchmark_stats          = BENCHMARK.out.stats
    ch_multiqc_files            = ch_multiqc_files.mix(ch_benchmark_stats)
    ch_versions                 = ch_versions.mix(BENCHMARK.out.versions.first())

    // Collate and save software versions
    ch_versions                 = ch_versions
                                | unique
                                | filter { yml -> yml }

    ch_versions_yml             = softwareVersionsToYAML(ch_versions)
                                | collectFile(
                                    storeDir: "${params.outdir}/pipeline_info",
                                    name: 'genepal_software_mqc_versions.yml',
                                    sort: true,
                                    newLine: true,
                                    cache: false
                                )

    // MODULE: MultiQC
    ch_multiqc_config           = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    ch_workflow_summary         = Channel.value( paramsSummaryMultiqc ( paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json") ) )
                                | collectFile(name: 'workflow_summary_mqc.yaml')

    ch_methods_description      = ch_versions_yml
                                | map { versions_yml ->
                                    methodsDescriptionText ( file("$projectDir/assets/methods_description_template.yml", checkIfExists: true), versions_yml )
                                }
                                | collectFile(name: 'methods_description_mqc.yaml', sort: true)

    ch_multiqc_extra_files      = Channel.empty()
                                | mix(ch_workflow_summary)
                                | mix(ch_versions_yml)
                                | mix(ch_methods_description)

    MULTIQC (
        ch_multiqc_files
        | map { meta, file -> file }
        | mix(ch_multiqc_extra_files)
        | collect,
        ch_multiqc_config.toList(),
        [],
        [],
        [],
        []
    )

    ch_versions                 = ch_versions.mix(MULTIQC.out.versions)

    // MODULE: GENEPALREPORT
    ch_pipeline_info            = ch_workflow_summary
                                | mix(ch_versions_yml)
                                | mix(ch_methods_description)


    GENEPALREPORT (
        ch_saved_marked_gff3        .map { meta, file -> file   }   .collect()              ,
        ch_orthofinder_statistics   .map { meta, dir -> dir     }   .collect()  .ifEmpty([]),
        ch_orthofinder_hogs         .map { meta, dir -> dir     }   .collect()  .ifEmpty([]),
        ch_busco_fasta_summary      .map { meta, file -> file   }   .collect()  .ifEmpty([]),
        ch_busco_gff_summary        .map { meta, file -> file   }   .collect()  .ifEmpty([]),
        ch_benchmark_stats          .map { meta, file -> file   }   .collect()  .ifEmpty([]),
        ch_pipeline_info                                            .collect()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
