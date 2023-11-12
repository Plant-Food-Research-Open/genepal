include { GUNZIP as GUNZIP_TARGET_ASSEMBLY      } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TE_LIBRARY           } from '../../modules/nf-core/gunzip'
include { FASTA_VALIDATE                        } from '../../modules/local/fasta_validate'
include { REPEATMASKER                          } from '../../modules/kherronism/repeatmasker'
include { STAR_GENOMEGENERATE                   } from '../../modules/nf-core/star/genomegenerate'

include { FASTA_EDTA                            } from '../../subworkflows/local/fasta_edta'

workflow PREPARE_ASSEMBLY {
    take:
    target_assembly     // channel: [ meta, fasta ]
    te_library          // channel: [ meta, fasta ]

    main:
    // MODULE: GUNZIP_TARGET_ASSEMBLY
    target_assembly
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { tech_target_assembly_branch }

    GUNZIP_TARGET_ASSEMBLY(
        tech_target_assembly_branch.gz
    )
    .gunzip
    | mix(
        tech_target_assembly_branch.rest
    )
    | set { ch_gunzip_target_assembly }

    // MODULE: FASTA_VALIDATE
    FASTA_VALIDATE(ch_gunzip_target_assembly)
    .valid_fasta
    | set { ch_validated_target_assembly }

    // MODULE: GUNZIP_TE_LIBRARY
    te_library
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_te_library_branch }

    GUNZIP_TE_LIBRARY(
        ch_te_library_branch.gz
    )
    .gunzip
    | mix(
        ch_te_library_branch.rest
    )
    | set { ch_gunzip_te_library }

    // SUBWORKFLOW: FASTA_EDTA
    ch_validated_target_assembly
    | join(
        ch_gunzip_te_library, remainder: true
    )
    | filter { meta, assembly, teLib ->
        teLib == null
    }
    | map { meta, assembly, teLib -> [meta, assembly] }
    | FASTA_EDTA
    
    // MODULE: REPEATMASKER
    ch_validated_target_assembly
    | join(
        FASTA_EDTA.out.te_lib_fasta.mix(ch_gunzip_te_library)
    )
    | set { ch_assembly_n_te_lib }

    REPEATMASKER(
        ch_assembly_n_te_lib.map { meta, assembly, teLib -> [meta, assembly] },
        ch_assembly_n_te_lib.map { meta, assembly, teLib -> teLib },
    )

    // MODULE: STAR_GENOMEGENERATE
    def star_ignore_sjdbgtf = true
    STAR_GENOMEGENERATE(
        ch_validated_target_assembly,
        ch_validated_target_assembly.map { meta, maskedFasta -> [meta, []] },
        star_ignore_sjdbgtf
    )
    .index
    | set { ch_assembly_index }

    Channel.empty()
    | mix(FASTA_VALIDATE.out.versions.first())
    | mix(GUNZIP_TE_LIBRARY.out.versions.first())
    | mix(FASTA_EDTA.out.versions)
    | mix(REPEATMASKER.out.versions.first())
    | mix(STAR_GENOMEGENERATE.out.versions.first())
    | mix(GUNZIP_TARGET_ASSEMBLY.out.versions.first())
    | set { ch_versions }
    
    emit:
    target_assemby          = ch_validated_target_assembly  // channel: [ meta, fasta ]
    masked_target_assembly  = REPEATMASKER.out.fasta_masked // channel: [ meta, fasta ]
    target_assemby_index    = ch_assembly_index             // channel: [ meta, star_index ]
    versions                = ch_versions                   // channel: [ versions.yml ]
}