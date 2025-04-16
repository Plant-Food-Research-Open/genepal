include { AGAT_SPMERGEANNOTATIONS                       } from '../../modules/nf-core/agat/spmergeannotations/main'
include { GT_GFF3                                       } from '../../modules/nf-core/gt/gff3/main'
include { AGAT_SPFILTERBYORFSIZE as FILTER_BY_ORF_SIZE  } from '../../modules/nf-core/agat/spfilterbyorfsize/main'
include { AGAT_CONVERTSPGXF2GXF                         } from '../../modules/nf-core/agat/convertspgxf2gxf/main'

workflow GFF_MERGE_CLEANUP {
    take:
    ch_braker_gff               // Channel: [ meta, gff ]
    ch_liftoff_gff              // Channel: [ meta, gff ]
    val_filter_by_aa_length     // val(null|Integer)
    val_append_genome_prefix    // val(true|false)

    main:
    ch_versions                 = Channel.empty()

    ch_gff_branch               = ch_braker_gff
                                | join(ch_liftoff_gff, remainder:true)
                                | branch { _meta, braker_gff, liftoff_gff ->
                                    both        : (     braker_gff      &&      liftoff_gff )
                                    braker_only : (     braker_gff      && ( !  liftoff_gff ) )
                                    liftoff_only: ( ( ! braker_gff )    &&      liftoff_gff )
                                }

    // MODULE: AGAT_SPMERGEANNOTATIONS
    AGAT_SPMERGEANNOTATIONS(
        ch_gff_branch.both.map { meta, bg, lg -> [ meta, [ bg, lg ] ] },
        []
    )

    ch_merged_gff               = AGAT_SPMERGEANNOTATIONS.out.gff
                                | mix ( ch_gff_branch.liftoff_only.map { meta, _braker_gff, liftoff_gff -> [ meta, liftoff_gff ] } )
                                | mix ( ch_gff_branch.braker_only.map { meta, braker_gff, _liftoff_gff -> [ meta, braker_gff ] } )
    ch_versions                 = ch_versions.mix(AGAT_SPMERGEANNOTATIONS.out.versions.first())

    // MODULE: AGAT_SPFILTERBYORFSIZE as FILTER_BY_ORF_SIZE
    ch_filter_input             = ch_merged_gff
                                | branch {
                                    filter: val_filter_by_aa_length != null
                                    pass: val_filter_by_aa_length == null
                                }

    FILTER_BY_ORF_SIZE ( ch_filter_input.filter, [] )

    ch_filtered_gff             = FILTER_BY_ORF_SIZE.out.passed_gff
                                | mix ( ch_filter_input.pass )
    ch_versions                 = ch_versions.mix(FILTER_BY_ORF_SIZE.out.versions.first())

    // MODULE: GT_GFF3
    GT_GFF3 ( ch_filtered_gff )

    ch_gt_gff                   = GT_GFF3.out.gt_gff3
    ch_versions                 = ch_versions.mix(GT_GFF3.out.versions.first())

    // COLLECTFILE: Format GT_GFF3 output
    ch_gt_formatted_gff         = ch_gt_gff
                                | map { meta, gff ->

                                    def lines = gff.readLines()
                                        .collect { line ->
                                            if ( line.startsWith('##') ) { return line }
                                            if ( line.startsWith('#') ) { return '' }

                                            def cols    = line.split('\t')
                                            def program = cols[1]
                                            def feat    = cols[2]
                                            def atts    = cols[8]

                                            def atts_r  = ''
                                            // Remove attributes and use AGAT_CONVERTSPGXF2GXF
                                            // to create attributes based on sequential layout

                                            def feat_r  = feat == 'transcript' ? 'mRNA' : feat
                                            // Use mRNA inplace of transcript

                                            if ( feat_r != 'mRNA' || program != 'Liftoff' ) {
                                                return ( cols[0..1] + [ feat_r ] + cols[3..7] + [ atts_r ] ).join('\t')
                                            }

                                            def tx_id   = ( atts =~ /ID=([^;]*)/ )[0][1]
                                            def matches = ( atts =~ /liftoffID=([^;]*)/ )

                                            def liftoffID = matches ? matches[0][1] : tx_id

                                            def atts_g  = "liftoffID=$liftoffID"

                                            return ( cols[0..1] + [ feat_r ] + cols[3..7] + [ atts_g ] ).join('\t')

                                        }.join('\n')

                                    [ "${meta.id}.bare.gff" ] + [ lines ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.bare', '') ], file ]
                                }

    // MODULE: AGAT_CONVERTSPGXF2GXF
    AGAT_CONVERTSPGXF2GXF ( ch_gt_formatted_gff )

    ch_agat_gff                 = AGAT_CONVERTSPGXF2GXF.out.output_gff
    ch_versions                 = ch_versions.mix(AGAT_CONVERTSPGXF2GXF.out.versions.first())

    // COLLECTFILE: Format AGAT_CONVERTSPGXF2GXF output and only allow: [ 'gene', 'mRNA', 'exon', 'CDS' ]
    ch_agat_formatted_gff       = ch_agat_gff
                                | map { meta, gff ->

                                    def filtered_lines = gff.readLines()
                                        .findAll { line ->
                                            if ( line.startsWith('#') ) { return true }

                                            def cols    = line.split('\t')
                                            def feat    = cols[2].trim()

                                            ( feat in [ 'gene', 'mRNA', 'exon', 'CDS' ] )
                                            ? true
                                            : false
                                        }
                                        .collect { line ->
                                            if ( line.startsWith('#') ) { return line }

                                            def cols    = line.split('\t')
                                            def program = cols[1]
                                            def feat    = cols[2]
                                            def atts    = cols[8]
                                            def atts_r  = atts.replace('-', '').replace('agat', '')

                                            if ( feat != 'mRNA' || program != 'Liftoff' ) {
                                                return ( cols[0..7] + [ atts_r ] ).join('\t')
                                            }

                                            def oldID   = ( atts =~ /liftoffID=([^;]*)/ )[0][1]
                                            def newID   = ( atts =~ /ID=([^;]*)/ )[0][1].replace('-', '').replace('agat', '')
                                            def pID     = ( atts =~ /Parent=([^;]*)/ )[0][1].replace('-', '').replace('agat', '')
                                            def atts_g  = "ID=${newID};Parent=${pID};liftoffID=${oldID}"

                                            return ( cols[0..7] + [ atts_g ] ).join('\t')
                                        }

                                    def tx_formatted_lines  = []
                                    def gene_counter        = 0
                                    def current_gene_id     = ''
                                    def current_mrna_id     = -1
                                    def current_exon_id     = -1
                                    def current_cds_id      = -1

                                    filtered_lines.each { line ->
                                        if (line.startsWith('#')) {
                                            tx_formatted_lines << line
                                            return
                                        }
                                        def cols    = line.split('\t')
                                        def feat    = cols[2]
                                        def atts    = cols[8]
                                        def id      = (atts =~ /ID=([^;]*)/)[0][1]

                                        if (feat == 'gene') {
                                            gene_counter += 1
                                            def new_gene_id = val_append_genome_prefix ? "${meta.id}.g${gene_counter}" : "g${gene_counter}"
                                            def updated_atts = atts.replaceAll(/ID=[^;]+/, "ID=${new_gene_id}")
                                            tx_formatted_lines << (cols[0..7] + [updated_atts]).join('\t')
                                            current_gene_id = new_gene_id
                                            current_mrna_id = 0
                                            return
                                        }

                                        if (feat == 'mRNA') {
                                            current_mrna_id += 1
                                            current_exon_id = 0
                                            current_cds_id = 0
                                            def matches = (atts =~ /liftoffID=([^;]*)/)
                                            def liftoffIDStr = matches ? ";liftoffID=${matches[0][1]}" : ''
                                            tx_formatted_lines << (
                                                (cols[0..7] + ["ID=${current_gene_id}.t${current_mrna_id};Parent=${current_gene_id}${liftoffIDStr}"]).join('\t')
                                            )
                                            return
                                        }

                                        if (feat == 'exon') {
                                            current_exon_id += 1
                                            tx_formatted_lines << (
                                                (cols[0..7] + ["ID=${current_gene_id}.t${current_mrna_id}.exon${current_exon_id};Parent=${current_gene_id}.t${current_mrna_id}"]).join('\t')
                                            )
                                            return
                                        }

                                        if (feat == 'CDS') {
                                            current_cds_id += 1
                                            tx_formatted_lines << (
                                                (cols[0..7] + ["ID=${current_gene_id}.t${current_mrna_id}.cds${current_cds_id};Parent=${current_gene_id}.t${current_mrna_id}"]).join('\t')
                                            )
                                            return
                                        }
                                    }

                                    [ "${meta.id}.agat.cleanup.gff" ] + [ tx_formatted_lines.join('\n') ]
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.agat.cleanup', '') ], file ]
                                }

    emit:
    gff                         = ch_agat_formatted_gff     // [ meta, gff ]
    versions                    = ch_versions               // [ versions.yml ]
}
