include { BRAKER3                               } from '../../modules/nf-core/braker3'
include { FILE_GUNZIP as BRAKER_GFF3_GUNZIP     } from '../../subworkflows/local/file_gunzip'
include { FILE_GUNZIP as BRAKER_HINTS_GUNZIP    } from '../../subworkflows/local/file_gunzip'
include { GFFREAD as FILTER_INVALID_ORFS        } from '../../modules/nf-core/gffread'

workflow FASTA_BRAKER3 {
    take:
    ch_valid_target_assembly    // channel: [ meta, fasta ]; meta ~ [ id: target_assembly ]; All input assemblies
    ch_masked_target_assembly   // channel: [ meta, fasta ]; Assemblies that don't have BRAKER annotations in the input sheet
    ch_rnaseq_bam               // channel: [ meta, bam ]
    ch_ext_prots_fasta          // channel: [ meta2, fasta ]; meta2 ~ [ id: ext_protein_seqs ]
    ch_braker_annotation        // channel: [ meta, gff3, hints.gff ]

    main:
    ch_versions                 = Channel.empty()


    ch_braker_inputs            = ch_masked_target_assembly
                                | join(ch_rnaseq_bam, remainder: true)
                                | combine(
                                    ch_ext_prots_fasta.map { meta, fasta -> fasta }.ifEmpty(null)
                                )
                                | map { meta, fasta, bam, prots -> [ meta, fasta, bam ?: [], prots ?: [] ] }

    def rnaseq_sets_dirs        = []
    def rnaseq_sets_ids         = []
    def hintsfile               = []

    // MODULE: BRAKER3
    BRAKER3(
        ch_braker_inputs.map { meta, fasta, bam, prots -> [meta, fasta] },
        ch_braker_inputs.map { meta, fasta, bam, prots -> bam },
        rnaseq_sets_dirs,
        rnaseq_sets_ids,
        ch_braker_inputs.map { meta, fasta, bam, prots -> prots },
        hintsfile
    )

    ch_braker_gff3              = BRAKER3.out.gff3
                                | mix( ch_braker_annotation.map { meta, gff3, hints -> [ meta, gff3 ] } )
    ch_braker_hints             = BRAKER3.out.hintsfile
                                | mix( ch_braker_annotation.map { meta, gff3, hints -> [ meta, hints ] } )
    ch_versions                 = ch_versions.mix(BRAKER3.out.versions.first())

    // WORKFLOW: FILE_GUNZIP as BRAKER_GFF3_GUNZIP
    BRAKER_GFF3_GUNZIP ( ch_braker_gff3 )

    ch_braker_gff3              = BRAKER_GFF3_GUNZIP.out.gunzip
    ch_versions                 = ch_versions.mix(BRAKER_GFF3_GUNZIP.out.versions)

    // WORKFLOW: FILE_GUNZIP as BRAKER_HINTS_GUNZIP
    BRAKER_HINTS_GUNZIP ( ch_braker_hints )

    ch_braker_hints             = BRAKER_HINTS_GUNZIP.out.gunzip
    ch_versions                 = ch_versions.mix(BRAKER_HINTS_GUNZIP.out.versions)


    // MODULE: GFFREAD as FILTER_INVALID_ORFS
    ch_filter_inputs            = ch_braker_gff3
                                | join(ch_valid_target_assembly)
                                | multiMap { meta, gff3, fasta ->
                                    gff:   [ meta, gff3 ]
                                    fasta: fasta
                                }

    FILTER_INVALID_ORFS (
        ch_filter_inputs.gff,
        ch_filter_inputs.fasta
    )

    ch_filtered_gff3            = FILTER_INVALID_ORFS.out.gffread_gff
    ch_versions                 = ch_versions.mix(FILTER_INVALID_ORFS.out.versions)

    emit:
    braker_gff3                 = ch_filtered_gff3      // [ meta, gff3 ]
    braker_hints                = ch_braker_hints       // [ meta, hints.gff ]
    versions                    = ch_versions           // [ versions.yml ]
}
