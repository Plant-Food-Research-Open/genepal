include { GUNZIP as GUNZIP_TARGET_ASSEMBLY      } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TE_LIBRARY           } from '../../modules/nf-core/gunzip'
include { FASTAVALIDATOR                        } from '../../modules/nf-core/fastavalidator'
include { REPEATMODELER_BUILDDATABASE           } from '../../modules/pfr/repeatmodeler/builddatabase'
include { REPEATMODELER_REPEATMODELER           } from '../../modules/pfr/repeatmodeler/repeatmodeler'
include { REPEATMASKER                          } from '../../modules/kherronism/repeatmasker'
include { STAR_GENOMEGENERATE                   } from '../../modules/nf-core/star/genomegenerate'

include { FASTA_EDTA_LAI                        } from '../../subworkflows/pfr/fasta_edta_lai'

workflow PREPARE_ASSEMBLY {
    take:
    target_assembly             // channel: [ meta, fasta ]
    te_library                  // channel: [ meta, fasta ]
    repeat_annotator            // val(String), 'repeatmodeler' or 'edta'
    exclude_assemblies          // channel: val(assembly_x,assembly_y)

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


    // MODULE: FASTAVALIDATOR
    FASTAVALIDATOR ( ch_gunzip_assembly )

    ch_validated_assembly       = ch_gunzip_assembly
                                | join(FASTAVALIDATOR.out.success_log)
                                | map { meta, fasta, log -> [ meta, fasta ] }
    ch_versions                 = ch_versions.mix(FASTAVALIDATOR.out.versions.first())

    FASTAVALIDATOR.out.error_log
    | map { meta, log ->
        System.err.println("WARNING: FASTAVALIDATOR failed for ${meta.id} with error: ${log}. ${meta.id} is excluded from further analysis.")
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
    ch_annotator_inputs         = ch_validated_assembly
                                | join(
                                    ch_gunzip_te_library, remainder: true
                                )
                                | filter { meta, assembly, teLib ->
                                    teLib == null
                                }
                                | map { meta, assembly, teLib -> [meta, assembly] }

    ch_edta_inputs              = repeat_annotator != 'edta'
                                ? Channel.empty()
                                : ch_annotator_inputs
                                | combine( exclude_assemblies )
                                | map { meta, fasta, ex_assemblies ->
                                    def ex_list = ex_assemblies.split(",")

                                    if ( !( ex_list.contains( meta.id ) ) ) {
                                        [ meta, fasta ]
                                    }
                                }

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
                                | combine( exclude_assemblies )
                                | map { meta, fasta, ex_assemblies ->
                                    def ex_list = ex_assemblies.split(",")

                                    if ( !( ex_list.contains( meta.id ) ) ) {
                                        [ meta, fasta ]
                                    }
                                }

    REPEATMODELER_BUILDDATABASE ( ch_repeatmodeler_inputs )

    ch_versions                 = ch_versions.mix(REPEATMODELER_BUILDDATABASE.out.versions.first())

    // MODULE: REPEATMODELER_REPEATMODELER
    REPEATMODELER_REPEATMODELER ( REPEATMODELER_BUILDDATABASE.out.db )

    ch_assembly_and_te_lib      = ch_validated_assembly
                                | join(
                                    repeat_annotator == 'edta'
                                    ? FASTA_EDTA_LAI.out.te_lib_fasta.mix(ch_gunzip_te_library)
                                    : REPEATMODELER_REPEATMODELER.out.fasta.mix(ch_gunzip_te_library)
                                )

    ch_versions                 = ch_versions.mix(REPEATMODELER_REPEATMODELER.out.versions.first())

    // MODULE: REPEATMASKER
    REPEATMASKER(
        ch_assembly_and_te_lib.map { meta, assembly, teLib -> [meta, assembly] },
        ch_assembly_and_te_lib.map { meta, assembly, teLib -> teLib },
    )

    ch_versions                 = ch_versions.mix(REPEATMASKER.out.versions.first())

    // MODULE: STAR_GENOMEGENERATE
    ch_genomegenerate_inputs    = ch_validated_assembly
                                | combine( exclude_assemblies )
                                | map { meta, fasta, ex_assemblies ->
                                    def ex_list = ex_assemblies.split(",")

                                    if ( !( ex_list.contains( meta.id ) ) ) {
                                        [ meta, fasta ]
                                    }
                                }


    STAR_GENOMEGENERATE(
        ch_genomegenerate_inputs,
        ch_genomegenerate_inputs.map { meta, fasta -> [ [], [] ] }
    )

    ch_assembly_index           = STAR_GENOMEGENERATE.out.index
    ch_versions                 = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

    emit:
    target_assemby              = ch_validated_assembly         // channel: [ meta, fasta ]
    masked_target_assembly      = REPEATMASKER.out.fasta_masked // channel: [ meta, fasta ]
    target_assemby_index        = ch_assembly_index             // channel: [ meta, star_index ]
    versions                    = ch_versions                   // channel: [ versions.yml ]
}
