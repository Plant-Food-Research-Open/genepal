include { AGAT_SPFILTERFEATUREFROMKILLLIST                  } from '../../modules/pfr/agat/spfilterfeaturefromkilllist/main'

workflow PURGE_NOHIT_MODELS {
    take:
    ch_target_gff               // [ meta, gff ]
    ch_eggnogmapper_hits        // [ meta, hits ]
    val_purge_nohits            // val(true|false)

    main:
    ch_versions                 = Channel.empty()

    // COLLECTFILE: Transcript level kill list
    ch_kill_list                = ch_target_gff
                                | join(ch_eggnogmapper_hits)
                                | map { meta, gff, hits ->

                                    def tx_with_hits = hits.readLines()
                                        .collect { it.split('\t')[0] }
                                        .sort(false)
                                        .unique()

                                    def tx_in_gff = gff.readLines()
                                        .findAll { line ->
                                            if ( line.startsWith('#') || line == '' ) { return false }

                                            def feat = line.split('\t')[2]
                                            ( feat == 'transcript' || feat == 'mRNA' )
                                        }
                                        .collect { it ->
                                            def attrs = it.split('\t')[8]

                                            ( attrs =~ /ID=([^;]*)/ )[0][1]
                                        }
                                        .sort(false)
                                        .unique()

                                    def tx_without_hits = tx_in_gff - tx_with_hits

                                    [ "${meta.id}.kill.list.txt" ] + tx_without_hits.join('\n')
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.kill.list', '') ], file ]
                                }

    // MODULE: AGAT_SPFILTERFEATUREFROMKILLLIST
    ch_agat_kill_inputs         = ! val_purge_nohits
                                ? Channel.empty()
                                : ch_target_gff
                                | join(ch_kill_list)


    AGAT_SPFILTERFEATUREFROMKILLLIST(
        ch_agat_kill_inputs.map { meta, gff, kill -> [ meta, gff ] },
        ch_agat_kill_inputs.map { meta, gff, kill -> kill },
        [] // default config
    )

    ch_target_purged_gff        = AGAT_SPFILTERFEATUREFROMKILLLIST.out.gff
    ch_versions                 = ch_versions.mix(AGAT_SPFILTERFEATUREFROMKILLLIST.out.versions.first())

    emit:
    purged_gff                  = ch_target_purged_gff.mix(val_purge_nohits ? Channel.empty() : ch_target_gff)
    versions                    = ch_versions   // [ versions.yml ]
}
