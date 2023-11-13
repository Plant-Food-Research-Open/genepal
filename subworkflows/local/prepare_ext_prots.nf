include { GUNZIP                        } from '../../modules/nf-core/gunzip'
include { CAT_CAT as CAT_PROTEIN_FASTAS } from '../../modules/nf-core/cat/cat'

workflow PREPARE_EXT_PROTS {
    take:
    ch_ext_prot_fastas          // Channel: [ meta, fasta ]
    
    main:
    ch_ext_prot_fastas
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_ext_prot_seqs_branch }

    // MODULE: GUNZIP
    GUNZIP(
        ch_ext_prot_seqs_branch.gz
    )
    .gunzip
    | mix(
        ch_ext_prot_seqs_branch.rest
    )
    | set { ch_ext_prot_gunzip_fastas }

    // MODULE: CAT_PROTEIN_FASTAS
    ch_ext_prot_gunzip_fastas
    | map { meta, filePath -> filePath }
    | collect
    | map { fileList -> [[id:"ext_protein_seqs"], fileList] }
    | CAT_PROTEIN_FASTAS

    Channel.empty()
    | mix(GUNZIP.out.versions.first())
    | mix(CAT_PROTEIN_FASTAS.out.versions)
    | set { ch_versions }
    
    emit:
    ext_prots_fasta = CAT_PROTEIN_FASTAS.out.file_out   // Channel: [ meta, fasta ]
    versions        = ch_versions                       // Channel: [ versions.yml ]
}