include { GUNZIP as GUNZIP_FASTA    } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF      } from '../../modules/nf-core/gunzip'
include { GFFREAD                   } from '../../modules/nf-core/gffread'
include { LIFTOFF                   } from '../../modules/pfr/liftoff'

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

    // MODULE: GFFREAD
    ch_gffread_inputs               = ch_xref_gunzip_gff
                                    | map { meta, gff ->
                                        [ gff.getSimpleName(), meta, gff ]
                                    } // For meta insertion later, remove when GFFREAD has meta

    GFFREAD ( ch_gffread_inputs.map { name, meta, gff -> gff } )

    ch_gffread_gff                  = GFFREAD.out.gffread_gff
                                    | map { gff -> [ gff.getSimpleName(), gff ] }
                                    | join(ch_gffread_inputs)
                                    | map { fid, gffread_gff, meta, gff -> [ meta, gffread_gff ] }
                                    // meta insertion

    ch_versions                     = ch_versions.mix(GFFREAD.out.versions.first())

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
        ch_liftoff_inputs.map { meta, target_fa, ref_fa, ref_gff -> ref_gff }
    )

    ch_liftoff_gff3                 = LIFTOFF.out.polished_gff3
                                    | map { meta, gff -> [ [ id: meta.target_assemby ], gff ] }
                                    | groupTuple

    ch_versions                     = ch_versions.mix(LIFTOFF.out.versions.first())

    emit:
    gff3        = ch_liftoff_gff3               // [ meta, [ gff3 ] ]
    versions    = ch_versions                   // [ versions.yml ]
}
