include { GFFREAD as GFF2FASTA_FOR_EGGNOGMAPPER } from '../../modules/nf-core/gffread/main'
include { EGGNOGMAPPER                          } from '../../modules/nf-core/eggnogmapper/main'

workflow GFF_EGGNOGMAPPER {
    take:
    ch_gff                      // Channel: [ meta, gff ]
    ch_fasta                    // Channel: [ meta, fasta ]
    db_folder                   // val(db_folder)

    main:
    // Versions
    ch_versions                 = Channel.empty()

    // MODULE: GFFREAD as GFF2FASTA_FOR_EGGNOGMAPPER
    ch_gffread_inputs           = ch_gff
                                | join(ch_fasta)

    GFF2FASTA_FOR_EGGNOGMAPPER(
        ch_gffread_inputs.map { meta, gff, _fasta -> [ meta, gff ] },
        ch_gffread_inputs.map { _meta, _gff, fasta -> fasta }
    )

    ch_gffread_fasta            = GFF2FASTA_FOR_EGGNOGMAPPER.out.gffread_fasta
    ch_versions                 = ch_versions.mix(GFF2FASTA_FOR_EGGNOGMAPPER.out.versions.first())


    ch_eggnogmapper_inputs      = ! db_folder
                                ? Channel.empty()
                                : ch_gffread_fasta
                                | combine(Channel.fromPath(db_folder))

    EGGNOGMAPPER(
        ch_eggnogmapper_inputs.map { meta, fasta, _db -> [ meta, fasta ] },
        [],
        ch_eggnogmapper_inputs.map { _meta, _fasta, db -> db },
        [ [], [] ]
    )

    ch_eggnogmapper_annotations = EGGNOGMAPPER.out.annotations
    ch_eggnogmapper_orthologs   = EGGNOGMAPPER.out.orthologs
    ch_eggnogmapper_hits        = EGGNOGMAPPER.out.hits
    ch_versions                 = ch_versions.mix(EGGNOGMAPPER.out.versions.first())

    emit:
    eggnogmapper_annotations    = ch_eggnogmapper_annotations
    eggnogmapper_orthologs      = ch_eggnogmapper_orthologs
    eggnogmapper_hits           = ch_eggnogmapper_hits
    versions                    = ch_versions
}
