include { AGAT_CONVERTSPGFF2GTF                             } from '../../modules/nf-core/agat/convertspgff2gtf/main'
include { TSEBRA                                            } from '../../modules/pfr/tsebra/main'
include { AGAT_CONVERTSPGXF2GXF                             } from '../../modules/nf-core/agat/convertspgxf2gxf/main'
include { GFFCOMPARE as COMPARE_BRAKER_TO_LIFTOFF           } from '../../modules/nf-core/gffcompare/main'
include { AGAT_SPFILTERFEATUREFROMKILLLIST                  } from '../../modules/pfr/agat/spfilterfeaturefromkilllist/main'
include { GFFCOMPARE as VALIDATE_PURGING_BY_AGAT            } from '../../modules/nf-core/gffcompare/main'
include { AGAT_SPMERGEANNOTATIONS as MERGE_BRAKER_LIFTOFF   } from '../../modules/pfr/agat/spmergeannotations/main'

workflow PURGE_BREAKER_MODELS {
    take:
    braker_gff3                 // [ meta, gff3 ]
    braker_hints                // [ meta, gff ]
    liftoff_gff3                // [ meta, gff3 ]
    tsebra_config               // val(tsebra_config)

    main:
    ch_versions                 = Channel.empty()

    // MODULE: AGAT_CONVERTSPGFF2GTF
    AGAT_CONVERTSPGFF2GTF ( braker_gff3 )

    ch_braker_gtf               = AGAT_CONVERTSPGFF2GTF.out.output_gtf
    ch_versions                 = ch_versions.mix(AGAT_CONVERTSPGFF2GTF.out.versions.first())

    // COLLECTFILE: Prepare for TSEBRA
    ch_tsebra_input_gtf         = ch_braker_gtf
                                | map { meta, gtf ->

                                    def lines = gtf.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('#') ) { return line }

                                            def cols = line.split('\t')
                                            def feat = cols[2]

                                            if ( ! ( feat in [ 'gene', 'transcript', 'mRNA' ] ) ) { return line }

                                            def atts    = cols[8]
                                            def matches = atts =~ /ID ([^;]*)/
                                            def id      = matches[0][1]

                                            def feat_format = ( feat == 'mRNA' ) ? 'transcript' : feat

                                            return ( cols[0..1] + [ feat_format ] + cols[3..7] + [ id ] ).join('\t')
                                        }.join('\n')

                                    [ "${meta.id}.clean.gtf" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace(".clean", "") ], file ]
                                }

    // MODULE: TSEBRA
    ch_tsebra_inputs            = ch_tsebra_input_gtf
                                | join(braker_hints)
                                | combine(Channel.fromPath(tsebra_config))
    TSEBRA(
        ch_tsebra_inputs.map { meta, gtf, gff, cfg -> [ meta, [ gtf ] ] },
        ch_tsebra_inputs.map { meta, gtf, gff, cfg -> [ gff ] },
        [],
        ch_tsebra_inputs.map { meta, gtf, gff, cfg -> cfg }
    )

    ch_tsebra_gtf               = TSEBRA.out.tsebra_gtf
    ch_versions                 = ch_versions.mix(TSEBRA.out.versions.first())

    // COLLECTFILE: Format TSEBRA output
    ch_tsebra_formatted_gtf     = ch_tsebra_gtf
                                | map { meta, gtf ->

                                    def lines = gtf.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('#') ) { return line }

                                            def cols    = line.split('\t')
                                            def atts_r  = ''
                                            // Remove attributes and use AGAT_CONVERTSPGXF2GXF
                                            // to create attributes based on sequential layout

                                            return ( cols[0..7] + [ atts_r ] ).join('\t')
                                        }.join('\n')

                                    [ "${meta.id}.gtf" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName ], file ]
                                }

    // MODULE: AGAT_CONVERTSPGXF2GXF
    AGAT_CONVERTSPGXF2GXF ( ch_tsebra_formatted_gtf )

    ch_tsebra_formatted_gff     = AGAT_CONVERTSPGXF2GXF.out.output_gff
    ch_versions                 = ch_versions.mix(AGAT_CONVERTSPGXF2GXF.out.versions.first())

    // COLLECTFILE: Format AGAT_CONVERTSPGXF2GXF output
    ch_tsebra_gff               = ch_tsebra_formatted_gff
                                | map { meta, gtf ->

                                    def lines = gtf.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('#') ) { return line }

                                            def cols    = line.split('\t')
                                            def atts_r  = cols[8].replaceAll('-', '').replaceAll('agat', '')

                                            return ( cols[0..7] + [ atts_r ] ).join('\t')
                                        }.join('\n')

                                    [ "${meta.id}.gff3" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName ], file ]
                                }

    // MODULE: GFFCOMPARE as COMPARE_BRAKER_TO_LIFTOFF
    ch_comparison_inputs        = ch_tsebra_gff
                                | join(liftoff_gff3)


    COMPARE_BRAKER_TO_LIFTOFF (
        ch_comparison_inputs.map { meta, braker, liftoff -> [ meta, braker ] },
        [ [], [], [] ],
        ch_comparison_inputs.map { meta, braker, liftoff -> [ meta, liftoff ] },
    )

    ch_tracking                 = COMPARE_BRAKER_TO_LIFTOFF.out.tracking
    ch_versions                 = ch_versions.mix(COMPARE_BRAKER_TO_LIFTOFF.out.versions.first())

    // COLLECTFILE: Transcript level kill list
    ch_kill_list                = ch_tracking
                                | map { meta, tracking ->

                                    def kept_lines = tracking.readLines()
                                        .findAll { line ->
                                            def cols = line.split('\t')

                                            ( cols[3] != 'u' ) && ( cols[3] != 'p' )
                                        }

                                    def tx_kill_list = kept_lines
                                        .collect { line ->
                                            def cols = line.split('\t')

                                            def matched = cols[4] =~ /q1:([^\|]+)\|([^\|]+)/

                                            matched[0][2].trim()
                                        }.join('\n')

                                    [ "${meta.id}.kill.list.txt" ] + tx_kill_list
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.kill.list', '') ], file ]
                                }

    // MODULE: AGAT_SPFILTERFEATUREFROMKILLLIST
    ch_agat_kill_inputs         = ch_tsebra_gff
                                | join(ch_kill_list)


    AGAT_SPFILTERFEATUREFROMKILLLIST(
        ch_agat_kill_inputs.map { meta, gff, kill -> [ meta, gff ] },
        ch_agat_kill_inputs.map { meta, gff, kill -> kill },
        [] // default config
    )

    ch_braker_purged_gff        = AGAT_SPFILTERFEATUREFROMKILLLIST.out.gff
    ch_versions                 = ch_versions.mix(AGAT_SPFILTERFEATUREFROMKILLLIST.out.versions.first())

    // Handle case where liftoff is not present
    ch_all_braker_gff           = ch_tsebra_gff
                                | join(ch_braker_purged_gff, remainder:true)
                                | map { meta, tsebra_gff, purged_gff ->
                                    if ( purged_gff ) { return [ meta, purged_gff ] }
                                    if ( tsebra_gff ) { return [ meta, tsebra_gff ] }
                                }

    emit:
    braker_purged_gff           = ch_all_braker_gff     // [ meta, gff3 ]
    versions                    = ch_versions           // [ versions.yml ]
}
