include { AGAT_SPMERGEANNOTATIONS               } from '../../modules/pfr/agat/spmergeannotations/main'
include { GT_GFF3                               } from '../../modules/nf-core/gt/gff3/main'
include { AGAT_CONVERTSPGXF2GXF                 } from '../../modules/nf-core/agat/convertspgxf2gxf/main'

workflow GFF_MERGE_CLEANUP {
    take:
    ch_braker_gff               // Channel: [ meta, gff ]
    ch_liftoff_gff              // Channel: [ meta, gff ]

    main:
    ch_versions                 = Channel.empty()

    ch_gff_branch               = ch_braker_gff
                                | join(ch_liftoff_gff, remainder:true)
                                | branch { meta, braker_gff, liftoff_gff ->
                                    both        : (     braker_gff      &&      liftoff_gff )
                                    braker_only : (     braker_gff      && ( !  liftoff_gff ) )
                                    liftoff_only: ( ( ! braker_gff )    &&      liftoff_gff )
                                }

    // MODULE: AGAT_SPMERGEANNOTATIONS
    AGAT_SPMERGEANNOTATIONS(
        ch_gff_branch.both.map { meta, bg, lg -> [ meta, [ bg, lg ] ] },
        []
    )

    ch_merged_gff               = AGAT_SPMERGEANNOTATIONS.out.gff
                                | mix(ch_gff_branch.liftoff_only)
                                | mix(ch_gff_branch.braker_only)
    ch_versions                 = ch_versions.mix(AGAT_SPMERGEANNOTATIONS.out.versions.first())

    // MODULE: GT_GFF3
    GT_GFF3 ( ch_merged_gff )

    ch_gt_gff                   = GT_GFF3.out.gt_gff3
    ch_versions                 = ch_versions.mix(GT_GFF3.out.versions.first())

    // COLLECTFILE: Format GT_GFF3 output
    ch_gt_formatted_gff         = ch_gt_gff
                                | map { meta, gff ->

                                    def lines = gff.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('##') ) { return line }
                                            if ( line.startsWith('#') ) { return '' }

                                            def cols    = line.split('\t')
                                            def program = cols[1]
                                            def feat    = cols[2]
                                            def atts    = cols[8]

                                            def atts_r  = ''
                                            // Remove attributes and use AGAT_CONVERTSPGXF2GXF
                                            // to create attributes based on sequential layout

                                            def feat_r  = feat == 'transcript' ? 'mRNA' : feat
                                            // Use mRNA inplace of transcript

                                            if ( feat != 'gene' || program != 'Liftoff' ) {
                                                return ( cols[0..1] + [ feat_r ] + cols[3..7] + [ atts_r ] ).join('\t')
                                            }

                                            def gene_id = ( atts =~ /ID=([^;]*)/ )[0][1]
                                            def atts_g  = "liftoffID=$gene_id"

                                            return ( cols[0..7] + [ atts_g ] ).join('\t')

                                        }.join('\n')

                                    [ "${meta.id}.bare.gff" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.bare', '') ], file ]
                                }

    // MODULE: AGAT_CONVERTSPGXF2GXF
    AGAT_CONVERTSPGXF2GXF ( ch_gt_formatted_gff )

    ch_agat_gff                 = AGAT_CONVERTSPGXF2GXF.out.output_gff
    ch_versions                 = ch_versions.mix(AGAT_CONVERTSPGXF2GXF.out.versions.first())

    // COLLECTFILE: Format AGAT_CONVERTSPGXF2GXF output
    ch_final_gff                = ch_agat_gff
                                | map { meta, gff ->

                                    def lines = gff.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('#') ) { return line }

                                            def cols    = line.split('\t')
                                            def program = cols[1]
                                            def feat    = cols[2]
                                            def atts    = cols[8]
                                            def atts_r  = atts.replace('-', '').replace('agat', '')

                                            if ( feat != 'gene' || program != 'Liftoff' ) {
                                                return ( cols[0..7] + [ atts_r ] ).join('\t')
                                            }

                                            def oldID   = ( atts =~ /liftoffID=([^;]*)/ )[0][1]
                                            def newID   = ( atts =~ /ID=([^;]*)/ )[0][1].replace('-', '').replace('agat', '')
                                            def atts_g  = "ID=${newID};liftoffID=${oldID}"

                                            return ( cols[0..7] + [ atts_g ] ).join('\t')
                                        }.join('\n')

                                    [ "${meta.id}.agat.cleanup.gff" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.agat.cleanup', '') ], file ]
                                }

    emit:
    gff                         = ch_final_gff      // [ meta, gff ]
    versions                    = ch_versions       // [ versions.yml ]
}
