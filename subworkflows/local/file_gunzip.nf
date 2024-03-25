include { GUNZIP            } from '../../modules/nf-core/gunzip'

workflow FILE_GUNZIP {
    take:
    ch_input                // channel [ meta, archive ]

    main:
    ch_input_branch         = ch_input
                            | branch { meta, archive ->
                                gz: "$archive".endsWith('.gz')
                                rest: ! "$archive".endsWith('.gz')
                            }

    // MODULE: GUNZIP
    GUNZIP ( ch_input_branch.gz )

    emit:
    versions                = GUNZIP.out.versions.first()
    gunzip                  = GUNZIP.out.gunzip.mix( ch_input_branch.rest )
}
