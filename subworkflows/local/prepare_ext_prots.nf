include { GUNZIP                        } from '../../modules/nf-core/gunzip'
include { CAT_CAT as CAT_PROTEIN_FASTAS } from '../../modules/nf-core/cat/cat'

workflow PREPARE_EXT_PROTS {
    take:
    ch_ext_prot_fastas          // Channel: [ meta, fasta ]

    main:
    ch_versions                 = Channel.empty()

    // MODULE: GUNZIP
    ch_ext_prot_seqs_branch     = ch_ext_prot_fastas
                                | branch { meta, file ->
                                    gz: "$file".endsWith(".gz")
                                    rest: !"$file".endsWith(".gz")
                                }

    GUNZIP ( ch_ext_prot_seqs_branch.gz )

    ch_ext_prot_gunzip_fastas   = GUNZIP.out.gunzip.mix(ch_ext_prot_seqs_branch.rest)
                                | map { meta, filePath -> filePath }
                                | collect
                                | map { fileList -> [ [ id: "ext_protein_seqs" ], fileList ] }

    ch_versions                 = ch_versions.mix(GUNZIP.out.versions.first())

    // MODULE: CAT_CAT as CAT_PROTEIN_FASTAS
    CAT_PROTEIN_FASTAS ( ch_ext_prot_gunzip_fastas )

    ch_versions                 = ch_versions.mix(CAT_PROTEIN_FASTAS.out.versions)

    emit:
    ext_prots_fasta             = CAT_PROTEIN_FASTAS.out.file_out   // Channel: [ meta, fasta ]
    versions                    = ch_versions                       // Channel: [ versions.yml ]
}
