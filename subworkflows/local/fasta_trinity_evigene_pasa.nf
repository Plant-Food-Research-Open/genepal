include { TRINITY               } from '../../modules/nf-core/trinity'
include { EVIGENE_TR2AACDS      } from '../../modules/pfr/evigene/tr2aacds'

workflow FASTA_TRINITY_EVIGENE_PASA {
    take:
    reads_target                // channel: [ meta, assembly_id ]
    trim_reads                  // channel: [ meta, [ fq ] ]
    assembly_fasta              // channel: [ meta2, fasta ]

    main:
    ch_versions                 = Channel.empty()

    // MODULE: TRINITY
    ch_all_trim_reads           = trim_reads
                                | map { meta, fq ->
                                    [ [ id: 'all_trim_reads', single_end: meta.single_end ], fq ]
                                }
                                | groupTuple
                                | map { meta, group -> [ meta, group.flatten() ] }
                                // TODO: This logic might break when single_end
                                // true and false are mixed

    TRINITY ( ch_all_trim_reads )

    ch_trinity_fasta            = TRINITY.out.transcript_fasta
    ch_versions                 = ch_versions.mix(TRINITY.out.versions.first())

    // MODULE: EVIGENE_TR2AACDS
    EVIGENE_TR2AACDS ( ch_trinity_fasta )

    ch_evigene_okayset          = EVIGENE_TR2AACDS.out.okayset
    ch_versions                 = ch_versions.mix(EVIGENE_TR2AACDS.out.versions.first())

    emit:
    proteins                    = ch_evigene_okayset.map { meta, set -> [ meta, [ set.find { "$it" ==~ /.*\.okay\.aa/ } ] ] }
    versions                    = ch_versions
}
