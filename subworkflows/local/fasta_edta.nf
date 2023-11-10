include { SHORTEN_EDTA_IDS  } from '../../modules/local/edta/shorten_edta_ids'
include { EDTA              } from '../../modules/local/edta/edta'
include { RESTORE_EDTA_IDS  } from '../../modules/local/edta/restore_edta_ids'

workflow FASTA_EDTA {
    take:
    genome_fasta    // channel: [ meta, /path/fasta ]
    
    main:
    SHORTEN_EDTA_IDS(genome_fasta)
    .renamed_ids_fasta
    | EDTA

    RESTORE_EDTA_IDS(
        EDTA.out.te_lib_fasta,
        EDTA.out.intact_gff3.map { it[1] },
        EDTA.out.pass_list.map { it[1] },
        EDTA.out.out_file.map { it[1] },
        EDTA.out.te_anno_gff3.map { it[1] },
        SHORTEN_EDTA_IDS.out.renamed_ids_tsv.map { it[1] }
    )

    Channel.empty()
    | mix(
        SHORTEN_EDTA_IDS.out.versions.first()
    )
    | mix(
        EDTA.out.versions.first()
    )
    | mix(
        RESTORE_EDTA_IDS.out.versions.first()
    )
    | set { ch_versions }
    
    emit:
    te_lib_fasta    = RESTORE_EDTA_IDS.out.te_lib_fasta     // channel: [ meta, /path/fasta ]
    intact_gff3     = RESTORE_EDTA_IDS.out.intact_gff3      // channel: [ meta, /path/gff3 ]
    pass_list       = RESTORE_EDTA_IDS.out.pass_list        // channel: [ meta, /path/pass.list ]
    out_file        = RESTORE_EDTA_IDS.out.out_file         // channel: [ meta, /path/out.file ]
    te_anno_gff3    = RESTORE_EDTA_IDS.out.te_anno_gff3     // channel: [ meta, /path/gff3 ]
    renamed_ids_tsv = RESTORE_EDTA_IDS.out.renamed_ids_tsv  // channel: [ meta, /path/tsv ]
    versions        = ch_versions                           // channel: [ versions.yml ]
}