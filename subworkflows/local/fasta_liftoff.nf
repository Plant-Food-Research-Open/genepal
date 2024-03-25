include { GUNZIP as GUNZIP_FASTA                    } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF                      } from '../../modules/nf-core/gunzip/main'
include { GFFREAD as GFFREAD_BEFORE_LIFTOFF         } from '../../modules/nf-core/gffread/main'
include { LIFTOFF                                   } from '../../modules/nf-core/liftoff/main'
include { GFFCOMPARE as COMBINE_LIFTOFF_ANNOTATIONS } from '../../modules/nf-core/gffcompare/main'
include { GFFREAD as GFFREAD_AFTER_LIFTOFF          } from '../../modules/nf-core/gffread/main'

workflow FASTA_LIFTOFF {
    take:
    target_assemby                  // Channel: [ meta, fasta ]
    xref_fasta                      // Channel: [ meta2, fasta ]
    xref_gff                        // Channel: [ meta2, gff3 ]

    main:
    ch_versions                     = Channel.empty()

    // MODULE: GUNZIP as GUNZIP_FASTA
    ch_xref_fasta_branch            = xref_fasta
                                    | branch { meta, file ->
                                        gz: "$file".endsWith(".gz")
                                        rest: !"$file".endsWith(".gz")
                                    }

    GUNZIP_FASTA ( ch_xref_fasta_branch.gz )

    ch_xref_gunzip_fasta            = GUNZIP_FASTA.out.gunzip
                                    | mix(
                                        ch_xref_fasta_branch.rest
                                    )

    ch_versions                     = ch_versions.mix(GUNZIP_FASTA.out.versions.first())

    // MODULE: GUNZIP as GUNZIP_GFF
    ch_xref_gff_branch              = xref_gff
                                    | branch { meta, file ->
                                        gz: "$file".endsWith(".gz")
                                        rest: !"$file".endsWith(".gz")
                                    }

    GUNZIP_GFF ( ch_xref_gff_branch.gz )

    ch_xref_gunzip_gff              = GUNZIP_GFF.out.gunzip
                                    | mix(
                                        ch_xref_gff_branch.rest
                                    )

    ch_versions                     = ch_versions.mix(GUNZIP_GFF.out.versions.first())

    // MODULE: GFFREAD as GFFREAD_BEFORE_LIFTOFF
    ch_gffread_inputs               = ch_xref_gunzip_gff
                                    | map { meta, gff ->
                                        [ gff.baseName, meta, gff ]
                                    } // For meta insertion later, remove when GFFREAD has meta

    GFFREAD_BEFORE_LIFTOFF ( ch_gffread_inputs.map { name, meta, gff -> gff } )

    ch_gffread_gff                  = GFFREAD_BEFORE_LIFTOFF.out.gffread_gff
                                    | map { gff -> [ gff.baseName - '.gffread', gff ] }
                                    | join(ch_gffread_inputs)
                                    | map { fid, gffread_gff, meta, gff -> [ meta, gffread_gff ] }
                                    // meta insertion

    ch_versions                     = ch_versions.mix(GFFREAD_BEFORE_LIFTOFF.out.versions.first())

    // MODULE: LIFTOFF
    ch_liftoff_inputs               = target_assemby
                                    | combine(
                                        ch_xref_gunzip_fasta
                                        | join(
                                            ch_gffread_gff
                                        )
                                    )
                                    | map { meta, target_fa, ref_meta, ref_fa, ref_gff ->
                                        [
                                            [
                                                id: "${meta.id}.from.${ref_meta.id}",
                                                target_assemby: meta.id
                                            ],
                                            target_fa,
                                            ref_fa,
                                            ref_gff
                                        ]
                                    }

    LIFTOFF(
        ch_liftoff_inputs.map { meta, target_fa, ref_fa, ref_gff -> [ meta, target_fa ] },
        ch_liftoff_inputs.map { meta, target_fa, ref_fa, ref_gff -> ref_fa },
        ch_liftoff_inputs.map { meta, target_fa, ref_fa, ref_gff -> ref_gff },
        []
    )

    ch_liftoff_gff3                 = LIFTOFF.out.polished_gff3
                                    | map { meta, gff -> [ [ id: meta.target_assemby ], gff ] }
                                    | groupTuple

    ch_versions                     = ch_versions.mix(LIFTOFF.out.versions.first())

    // MODULE: GFFCOMPARE as COMBINE_LIFTOFF_ANNOTATIONS
    ch_gffcompare_inputs            = ch_liftoff_gff3
                                    | branch { meta, list_polished ->
                                        one: list_polished.size() == 1
                                        many: list_polished.size() > 1
                                    }

    COMBINE_LIFTOFF_ANNOTATIONS(
        ch_gffcompare_inputs.many,
        [ [], [], [] ],
        [ [], [] ],
    )

    ch_combined_gtf                 = COMBINE_LIFTOFF_ANNOTATIONS.out.combined_gtf
    ch_versions                     = ch_versions.mix(COMBINE_LIFTOFF_ANNOTATIONS.out.versions.first())

    // MODULE: GFFREAD as GFFREAD_AFTER_LIFTOFF
    ch_post_gffread_inputs          = ch_combined_gtf
                                    | map { meta, gtf ->
                                        [ gtf.baseName, meta, gtf ]
                                    } // For meta insertion later, remove when GFFREAD has meta

    GFFREAD_AFTER_LIFTOFF ( ch_post_gffread_inputs.map { name, meta, gtf -> gtf } )

    ch_combined_gff3                = GFFREAD_AFTER_LIFTOFF.out.gffread_gff
                                    | map { gff -> [ gff.baseName - '.gffread', gff ] }
                                    | join(ch_post_gffread_inputs)
                                    | map { fid, gffread_gff, meta, gtf -> [ meta, gffread_gff ] }
                                    // meta insertion
                                    | mix(ch_gffcompare_inputs.one)

    ch_versions                     = ch_versions.mix(GFFREAD_AFTER_LIFTOFF.out.versions.first())

    emit:
    gff3        = ch_combined_gff3              // [ meta, gff3 ]
    versions    = ch_versions                   // [ versions.yml ]
}
