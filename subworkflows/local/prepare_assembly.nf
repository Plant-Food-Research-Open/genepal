include { GUNZIP as GUNZIP_TARGET_ASSEMBLY      } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TE_LIBRARY           } from '../../modules/nf-core/gunzip'
include { SEQKIT as SEQKIT_FILTER               } from '../../modules/nf-core/seqkit/main.nf'
include { SEQKIT_RMDUP                          } from '../../modules/nf-core/seqkit/rmdup/main.nf'
include { FASTAVALIDATOR                        } from '../../modules/nf-core/fastavalidator'
include { REPEATMODELER_BUILDDATABASE           } from '../../modules/nf-core/repeatmodeler/builddatabase'
include { REPEATMODELER_REPEATMODELER           } from '../../modules/nf-core/repeatmodeler/repeatmodeler'
include { REPEATMASKER_REPEATMASKER             } from '../../modules/gallvp/repeatmasker/repeatmasker'
include { CUSTOM_RMOUTTOGFF3                    } from '../../modules/gallvp/custom/rmouttogff3'
include { STAR_GENOMEGENERATE                   } from '../../modules/nf-core/star/genomegenerate'

include { FASTA_EDTA_LAI                        } from '../../subworkflows/gallvp/fasta_edta_lai'

workflow PREPARE_ASSEMBLY {
    take:
    target_assembly             // channel: [ meta, fasta ]
    te_library                  // channel: [ meta, fasta ]
    repeat_annotator            // val(String), 'repeatmodeler' or 'edta'
    repeatmasker_save_outputs   // val(true/false)
    exclude_assemblies          // channel: val(assembly_x,assembly_y)
    ch_is_masked                // channel: [ meta, val(true|false) ]

    main:
    ch_versions                 = Channel.empty()

    // MODULE: GUNZIP_TARGET_ASSEMBLY
    target_assembly_branch      = target_assembly
                                | branch { meta, file ->
                                    gz: "$file".endsWith(".gz")
                                    rest: !"$file".endsWith(".gz")
                                }

    GUNZIP_TARGET_ASSEMBLY ( target_assembly_branch.gz )

    ch_gunzip_assembly          = GUNZIP_TARGET_ASSEMBLY.out.gunzip
                                | mix(
                                    target_assembly_branch.rest
                                )
    ch_versions                 = ch_versions.mix(GUNZIP_TARGET_ASSEMBLY.out.versions.first())

    // MODULE: SEQKIT_FILTER
    SEQKIT_FILTER ( ch_gunzip_assembly, params.min_contig_length )

    ch_filtered_assembly = SEQKIT_FILTER.out.fastx

    // MODULE: SEQKIT_RMDUP
    SEQKIT_RMDUP ( ch_filtered_assembly )

    ch_nondup_fw_assembly       = SEQKIT_RMDUP.out.log
                                | join(SEQKIT_RMDUP.out.fastx)
                                | map { meta, error_log, fasta ->
                                    if ( error_log.text.contains('0 duplicated records removed') ) {
                                        return [ meta, fasta ]
                                    }

                                    log.warn "FASTA validation failed for ${meta.id} due to presence of duplicate sequences.\n" +
                                        "${meta.id} is excluded from further analysis."

                                    return null
                                } // Fixed width assembly fasta without duplicates

    ch_versions                 = ch_versions.mix(SEQKIT_RMDUP.out.versions.first())

    // MODULE: FASTAVALIDATOR
    FASTAVALIDATOR ( ch_nondup_fw_assembly )

    ch_validated_assembly       = ch_nondup_fw_assembly
                                | join(FASTAVALIDATOR.out.success_log)
                                | map { meta, fasta, log -> [ meta, fasta ] }
    ch_versions                 = ch_versions.mix(FASTAVALIDATOR.out.versions.first())

    FASTAVALIDATOR.out.error_log
    | map { meta, log ->
        log.warn "FASTAVALIDATOR failed for ${meta.id} with error: ${log}. ${meta.id} is excluded from further analysis."
    }

    // MODULE: GUNZIP_TE_LIBRARY
    ch_te_library_branch        = te_library
                                | branch { meta, file ->
                                    gz: "$file".endsWith(".gz")
                                    rest: !"$file".endsWith(".gz")
                                }

    GUNZIP_TE_LIBRARY ( ch_te_library_branch.gz )

    ch_gunzip_te_library        = GUNZIP_TE_LIBRARY.out.gunzip
                                | mix(
                                    ch_te_library_branch.rest
                                )
    ch_versions                 = ch_versions.mix(GUNZIP_TE_LIBRARY.out.versions.first())

    // SUBWORKFLOW: FASTA_EDTA_LAI
    ch_unmasked_masked_branch   = ch_validated_assembly
                                | combine( exclude_assemblies )
                                | map { meta, fasta, ex_assemblies ->
                                    ex_assemblies.tokenize(",").contains( meta.id )
                                    ? null
                                    : [ meta, fasta ]
                                }
                                | join(
                                    ch_is_masked
                                )
                                | branch { meta, fasta, is_masked ->
                                    unmasked: ! is_masked
                                        return [ meta, fasta ]
                                    masked: is_masked
                                        return [ meta, fasta ]
                                }

    ch_annotator_inputs         = ch_unmasked_masked_branch.unmasked
                                | join(
                                    ch_gunzip_te_library, remainder: true
                                )
                                | filter { meta, assembly, teLib ->
                                    teLib == null && ( assembly != null )
                                }
                                | map { meta, assembly, teLib -> [ meta, assembly ] }

    ch_edta_inputs              = repeat_annotator != 'edta'
                                ? Channel.empty()
                                : ch_annotator_inputs

    FASTA_EDTA_LAI(
        ch_edta_inputs,
        [],
        true // Skip LAI
    )

    ch_versions                 = ch_versions.mix(FASTA_EDTA_LAI.out.versions.first())

    // MODULE: REPEATMODELER_BUILDDATABASE
    ch_repeatmodeler_inputs     = repeat_annotator != 'repeatmodeler'
                                ? Channel.empty()
                                : ch_annotator_inputs

    REPEATMODELER_BUILDDATABASE ( ch_repeatmodeler_inputs )

    ch_versions                 = ch_versions.mix(REPEATMODELER_BUILDDATABASE.out.versions.first())

    // MODULE: REPEATMODELER_REPEATMODELER
    REPEATMODELER_REPEATMODELER ( REPEATMODELER_BUILDDATABASE.out.db )

    ch_assembly_and_te_lib      = ch_unmasked_masked_branch.unmasked
                                | join(
                                    repeat_annotator == 'edta'
                                    ? FASTA_EDTA_LAI.out.te_lib_fasta.mix(ch_gunzip_te_library)
                                    : REPEATMODELER_REPEATMODELER.out.fasta.mix(ch_gunzip_te_library)
                                )

    ch_versions                 = ch_versions.mix(REPEATMODELER_REPEATMODELER.out.versions.first())

    // MODULE: REPEATMASKER_REPEATMASKER
    REPEATMASKER_REPEATMASKER(
        ch_assembly_and_te_lib.map { meta, assembly, teLib -> [ meta, assembly ] },
        ch_assembly_and_te_lib.map { meta, assembly, teLib -> teLib },
    )

    ch_masked_assembly          = ch_unmasked_masked_branch.masked
                                | mix(REPEATMASKER_REPEATMASKER.out.masked)

    ch_repeatmasker_out         = REPEATMASKER_REPEATMASKER.out.out
    ch_versions                 = ch_versions.mix(REPEATMASKER_REPEATMASKER.out.versions.first())

    // MODULE: CUSTOM_RMOUTTOGFF3
    ch_RMOUTTOGFF3_input        = repeatmasker_save_outputs
                                ? ch_repeatmasker_out
                                : Channel.empty()
    CUSTOM_RMOUTTOGFF3 ( ch_RMOUTTOGFF3_input )

    ch_versions                 = ch_versions.mix(CUSTOM_RMOUTTOGFF3.out.versions.first())

    // MODULE: STAR_GENOMEGENERATE
    ch_genomegenerate_inputs    = ch_validated_assembly
                                | combine( exclude_assemblies )
                                | map { meta, fasta, ex_assemblies ->
                                    ex_assemblies.tokenize(",").contains( meta.id )
                                    ? null
                                    : [ meta, fasta ]
                                }


    STAR_GENOMEGENERATE(
        ch_genomegenerate_inputs,
        ch_genomegenerate_inputs.map { meta, fasta -> [ [], [] ] }
    )

    ch_assembly_index           = STAR_GENOMEGENERATE.out.index
    ch_versions                 = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

    emit:
    target_assemby              = ch_validated_assembly         // channel: [ meta, fasta ]
    masked_target_assembly      = ch_masked_assembly            // channel: [ meta, fasta ]
    target_assemby_index        = ch_assembly_index             // channel: [ meta, star_index ]
    versions                    = ch_versions                   // channel: [ versions.yml ]
}
