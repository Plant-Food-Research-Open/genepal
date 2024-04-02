include { GFFCOMPARE as COMPARE_BRAKER_TO_LIFTOFF           } from '../../modules/nf-core/gffcompare/main'
include { AGAT_SPFILTERFEATUREFROMKILLLIST                  } from '../../modules/pfr/agat/spfilterfeaturefromkilllist/main'
include { GFFCOMPARE as VALIDATE_PURGING_BY_AGAT            } from '../../modules/nf-core/gffcompare/main'
include { AGAT_SPMERGEANNOTATIONS as MERGE_BRAKER_LIFTOFF   } from '../../modules/pfr/agat/spmergeannotations/main'

workflow MERGE_ANNOTATIONS {
    take:
    braker_gff3                 // [ meta, gff3 ]
    liftoff_gff3                // [ meta, gff3 ]

    main:
    ch_versions                 = Channel.empty()

    // MODULE: GFFCOMPARE as COMPARE_BRAKER_TO_LIFTOFF
    ch_comparison_inputs        = braker_gff3
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
                                        }

                                    [ "${meta.id}.kill.list.txt" ] + tx_kill_list
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace(".kill.list", "") ], file ]
                                }

    // MODULE: AGAT_SPFILTERFEATUREFROMKILLLIST
    ch_agat_kill_inputs         = braker_gff3
                                | join(ch_kill_list)


    AGAT_SPFILTERFEATUREFROMKILLLIST(
        ch_agat_kill_inputs.map { meta, gff, kill -> [ meta, gff ] },
        ch_agat_kill_inputs.map { meta, gff, kill -> kill },
        [] // default config
    )

    ch_braker_purged            = AGAT_SPFILTERFEATUREFROMKILLLIST.out.gff
    ch_versions                 = ch_versions.mix(AGAT_SPFILTERFEATUREFROMKILLLIST.out.versions.first())

    // MODULE: GFFCOMPARE as VALIDATE_PURGING_BY_AGAT
    ch_validation_inputs        = ch_braker_purged
                                | join(liftoff_gff3)

    VALIDATE_PURGING_BY_AGAT(
        ch_validation_inputs.map { meta, braker, liftoff -> [ meta, braker ] },
        [ [], [], [] ],
        ch_validation_inputs.map { meta, braker, liftoff -> [ meta, liftoff ] }
    )

    ch_validation_tracking      = VALIDATE_PURGING_BY_AGAT.out.tracking
    ch_versions                 = ch_versions.mix(VALIDATE_PURGING_BY_AGAT.out.versions.first())

    // VALIDATE: All gffcompare codes in [ 'u', 'p' ]
    ch_validation_output        = ch_validation_tracking
                                | map { meta, tracking ->

                                    def tx_list = tracking.readLines()
                                        .collect { line ->
                                            def cols    = line.split('\t')
                                            def matched = cols[4] =~ /q1:([^\|]+)\|([^\|]+)/
                                            def tx_id   = matched[0][2].trim()
                                            def code    = cols[3]

                                            if ( code != 'u' && code != 'p' ) {
                                                log.warn "In ${meta.id}, transcript $tx_id has code ${code}." +
                                                " This issue can break the downstream logic in the pipeline." +
                                                " Please contact the devs!"
                                            }

                                            [ tx_id, code ]
                                        }

                                    [ meta, tx_list ]
                                }

    // MODULE: AGAT_SPMERGEANNOTATIONS as MERGE_BRAKER_LIFTOFF
    ch_merge_inputs             = ch_braker_purged
                                | join(liftoff_gff3)
    MERGE_BRAKER_LIFTOFF(
        ch_merge_inputs.map { meta, braker, liftoff -> [ meta, [ braker, liftoff ] ] },
        []
    )

    ch_merged_gff               = MERGE_BRAKER_LIFTOFF.out.gff
    ch_versions                 = ch_versions.mix(MERGE_BRAKER_LIFTOFF.out.versions.first())

    emit:
    merged_gff                  = ch_merged_gff
    versions                    = ch_versions
}
