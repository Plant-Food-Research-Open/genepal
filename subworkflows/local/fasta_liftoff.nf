include { GUNZIP as GUNZIP_FASTA                                } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF                                  } from '../../modules/nf-core/gunzip/main'
include { GFFREAD as GFFREAD_BEFORE_LIFTOFF                     } from '../../modules/nf-core/gffread/main'
include { LIFTOFF                                               } from '../../modules/nf-core/liftoff/main'
include { AGAT_SPMERGEANNOTATIONS as MERGE_LIFTOFF_ANNOTATIONS  } from '../../modules/nf-core/agat/spmergeannotations/main'
include { AGAT_SPFLAGSHORTINTRONS                               } from '../../modules/gallvp/agat/spflagshortintrons/main'
include { AGAT_SPFILTERFEATUREFROMKILLLIST                      } from '../../modules/nf-core/agat/spfilterfeaturefromkilllist/main'
include { GFFREAD as GFFREAD_AFTER_LIFTOFF                      } from '../../modules/nf-core/gffread/main'
include { GT_GFF3 as GT_GFF3_AFTER_LIFTOFF                      } from '../../modules/nf-core/gt/gff3/main'
include { GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST                } from '../../subworkflows/local/gff_tsebra_spfilterfeaturefromkilllist'

workflow FASTA_LIFTOFF {
    take:
    target_assemby                  // Channel: [ meta, fasta ]
    xref_fasta                      // Channel: [ meta2, fasta ]
    xref_gff                        // Channel: [ meta2, gff3 ]
    val_filter_liftoff_by_hints     // val(true|false)
    braker_hints                    // [ meta, gff ]
    tsebra_config                   // Channel: [ cfg ]
    allow_isoforms                  // val(true|false)


    main:
    ch_versions                     = Channel.empty()

    // MODULE: GUNZIP as GUNZIP_FASTA
    ch_xref_fasta_branch            = xref_fasta
                                    | branch { meta, file ->
                                        gz: "$file".endsWith(".gz")
                                        rest: !"$file".endsWith(".gz")
                                    }

    GUNZIP_FASTA ( ch_xref_fasta_branch.gz )

    ch_xref_gunzip_fasta            = GUNZIP_FASTA.out.gunzip
                                    | mix(
                                        ch_xref_fasta_branch.rest
                                    )

    ch_versions                     = ch_versions.mix(GUNZIP_FASTA.out.versions.first())

    // MODULE: GUNZIP as GUNZIP_GFF
    ch_xref_gff_branch              = xref_gff
                                    | branch { meta, file ->
                                        gz: "$file".endsWith(".gz")
                                        rest: !"$file".endsWith(".gz")
                                    }

    GUNZIP_GFF ( ch_xref_gff_branch.gz )

    ch_xref_gunzip_gff              = GUNZIP_GFF.out.gunzip
                                    | mix(
                                        ch_xref_gff_branch.rest
                                    )

    ch_versions                     = ch_versions.mix(GUNZIP_GFF.out.versions.first())

    // MODULE: GFFREAD as GFFREAD_BEFORE_LIFTOFF
    GFFREAD_BEFORE_LIFTOFF ( ch_xref_gunzip_gff, [] )

    ch_gffread_gff                  = GFFREAD_BEFORE_LIFTOFF.out.gffread_gff
    ch_versions                     = ch_versions.mix(GFFREAD_BEFORE_LIFTOFF.out.versions.first())

    // MODULE: LIFTOFF
    ch_liftoff_inputs               = target_assemby
                                    | combine(
                                        ch_xref_gunzip_fasta
                                        | join(
                                            ch_gffread_gff
                                        )
                                    )
                                    | map { meta, target_fa, ref_meta, ref_fa, ref_gff ->
                                        [
                                            [
                                                id: "${meta.id}.from.${ref_meta.id}",
                                                target_assemby: meta.id
                                            ],
                                            target_fa,
                                            ref_fa,
                                            ref_gff
                                        ]
                                    }

    LIFTOFF(
        ch_liftoff_inputs.map { meta, target_fa, ref_fa, ref_gff -> [ meta, target_fa ] },
        ch_liftoff_inputs.map { meta, target_fa, ref_fa, ref_gff -> ref_fa },
        ch_liftoff_inputs.map { meta, target_fa, ref_fa, ref_gff -> ref_gff },
        []
    )

    ch_liftoff_gff3                 = LIFTOFF.out.polished_gff3
                                    | map { meta, gff -> [ [ id: meta.target_assemby ], gff ] }
                                    | groupTuple

    ch_versions                     = ch_versions.mix(LIFTOFF.out.versions.first())

    // MODULE: AGAT_SPMERGEANNOTATIONS as MERGE_LIFTOFF_ANNOTATIONS
    ch_merge_inputs                 = ch_liftoff_gff3
                                    | branch { meta, list_polished ->
                                        one: list_polished.size() == 1
                                        many: list_polished.size() > 1
                                    }

    MERGE_LIFTOFF_ANNOTATIONS(
        ch_merge_inputs.many,
        []
    )

    ch_merged_gff                   = MERGE_LIFTOFF_ANNOTATIONS.out.gff
                                    | mix(
                                        ch_merge_inputs.one
                                        | map { meta, gffs -> [ meta, gffs[0] ] }
                                        // Unlist the upstream groupTuple
                                    )
    ch_versions                     = ch_versions.mix(MERGE_LIFTOFF_ANNOTATIONS.out.versions.first())

    // MODULE: AGAT_SPFLAGSHORTINTRONS
    AGAT_SPFLAGSHORTINTRONS ( ch_merged_gff, [] )

    ch_flagged_gff                  = AGAT_SPFLAGSHORTINTRONS.out.gff
    ch_versions                     = ch_versions.mix(AGAT_SPFLAGSHORTINTRONS.out.versions.first())

    // COLLECTFILE: Kill list for valid_ORF=False transcripts
    // tRNA, rRNA, gene with any intron marked as
    // 'pseudo=' by AGAT/SPFLAGSHORTINTRONS
    ch_kill_list                    = ch_flagged_gff
                                    | map { meta, gff ->

                                        def tx_from_gff = gff.readLines()
                                            .findAll { it ->
                                                if ( it.startsWith('#') ) { return false }

                                                def cols = it.split('\t')
                                                def feat = cols[2]

                                                if ( feat in [ 'tRNA', 'rRNA' ] ) { return true }
                                                if ( feat !in [ 'transcript', 'mRNA', 'gene' ] ) { return false }

                                                def attrs = cols[8]
                                                ( attrs.contains('valid_ORF=False') || attrs.contains('pseudo=') )
                                            }
                                            .collect {
                                                def cols    = it.split('\t')
                                                def attrs   = cols[8]

                                                def matches = attrs =~ /ID=([^;]*)/

                                                return matches[0][1]
                                            }

                                        [ "${meta.id}.kill.list.txt" ] + tx_from_gff.join('\n')
                                    }
                                    | collectFile(newLine: true)
                                    | map { file ->
                                        [ [ id: file.baseName.replace('.kill.list', '') ], file ]
                                    }

    // MODULE: AGAT_SPFILTERFEATUREFROMKILLLIST
    ch_agat_kill_inputs             = ch_flagged_gff
                                    | join(ch_kill_list)


    AGAT_SPFILTERFEATUREFROMKILLLIST(
        ch_agat_kill_inputs.map { meta, gff, kill -> [ meta, gff ] },
        ch_agat_kill_inputs.map { meta, gff, kill -> kill },
        [] // default config
    )

    ch_liftoff_purged_gff           = AGAT_SPFILTERFEATUREFROMKILLLIST.out.gff
    ch_versions                     = ch_versions.mix(AGAT_SPFILTERFEATUREFROMKILLLIST.out.versions.first())

    // MODULE: GFFREAD as GFFREAD_AFTER_LIFTOFF
    GFFREAD_AFTER_LIFTOFF ( ch_liftoff_purged_gff, [] )

    ch_versions                     = ch_versions.mix(GFFREAD_AFTER_LIFTOFF.out.versions.first())

    // MODULE: GT_GFF3 as GT_GFF3_AFTER_LIFTOFF
    // To fix the phase mislabelled by LIFTOFF
    // https://github.com/Plant-Food-Research-Open/genepal/issues/115
    GT_GFF3_AFTER_LIFTOFF ( GFFREAD_AFTER_LIFTOFF.out.gffread_gff )

    ch_attr_trimmed_gff             = GT_GFF3_AFTER_LIFTOFF.out.gt_gff3
    ch_versions                     = ch_versions.mix(GT_GFF3_AFTER_LIFTOFF.out.versions.first())

    // SUBWORKFLOW: GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST
    GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST(
        val_filter_liftoff_by_hints ? ch_attr_trimmed_gff : Channel.empty(),
        braker_hints,
        tsebra_config,
        allow_isoforms,
        'liftoff'
    )

    ch_tsebra_killed_gff            = GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST.out.tsebra_killed_gff
    ch_versions                     = ch_versions.mix(GFF_TSEBRA_SPFILTERFEATUREFROMKILLLIST.out.versions)

    // Prepare output channel
    ch_output_gff                   = val_filter_liftoff_by_hints
                                    ? ch_tsebra_killed_gff
                                    : ch_attr_trimmed_gff

    emit:
    gff3                            = ch_output_gff         // [ meta, gff3 ]
    versions                        = ch_versions           // [ versions.yml ]
}
