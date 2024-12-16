include { GT_GFF3 as FINAL_GFF_CHECK    } from '../../modules/nf-core/gt/gff3/main'
include { GFFREAD as EXTRACT_PROTEINS   } from '../../modules/nf-core/gffread/main'
include { GFFREAD as EXTRACT_CDS        } from '../../modules/nf-core/gffread/main'
include { GFFREAD as EXTRACT_CDNA       } from '../../modules/nf-core/gffread/main'

workflow GFF_STORE {
    take:
    ch_target_gff               // [ meta, gff ]
    ch_eggnogmapper_annotations // [ meta, annotations ]
    ch_fasta                    // [ meta, fasta ]
    val_describe_gff            // val(true|false)

    main:
    ch_versions                 = Channel.empty()

    // COLLECTFILE: Add eggnogmapper hits to gff
    ch_described_gff            = ! val_describe_gff
                                ? Channel.empty()
                                : ch_target_gff
                                | join(ch_eggnogmapper_annotations)
                                | map { meta, gff, annotations ->
                                    def tx_annotations  = annotations.readLines()
                                        .findAll { ! it.startsWith('#') }
                                        .collect { line ->
                                            def cols    = line.split('\t')
                                            def id      = cols[0]
                                            def txt     = cols[7]
                                            def pfams   = cols[20]

                                            [ id, txt, pfams ]
                                        }
                                        .collect { id, txt, pfams ->
                                            if ( txt != '-' ) { return [ id, txt ] }
                                            if ( pfams != '-' ) { return [ id, "PFAMs: $pfams" ] }

                                            [ id, 'No eggnog description and PFAMs' ]
                                        }
                                        .collectEntries { id, txt ->
                                            [ id, txt ]
                                        }

                                    def gene_tx_annotations = [:]
                                    gff.readLines()
                                        .findAll { line ->
                                            if ( line.startsWith('#') || line == '' ) { return false }

                                            def cols    = line.split('\t')
                                            def feat    = cols[2]

                                            if ( ! ( feat == 'transcript' || feat == 'mRNA' ) ) { return false }

                                            return true
                                        }
                                        .each { line ->
                                            def cols    = line.split('\t')
                                            def atts    = cols[8]

                                            def matches = atts =~ /ID=([^;]*)/
                                            def tx_id   = matches[0][1]

                                            def matches_p= atts =~ /Parent=([^;]*)/
                                            def gene_id = matches_p[0][1]

                                            if ( ! gene_tx_annotations.containsKey(gene_id) ) {
                                                gene_tx_annotations[gene_id] = [:]
                                            }

                                            def anno    = tx_annotations.containsKey(tx_id)
                                                        ? java.net.URLEncoder.encode(tx_annotations[tx_id], "UTF-8").replace('+', '%20')
                                                        : java.net.URLEncoder.encode('Hypothetical protein | no eggnog hit', "UTF-8").replace('+', '%20')

                                            gene_tx_annotations[gene_id] += [ ( tx_id ): anno ]
                                        }

                                    gene_tx_annotations = gene_tx_annotations
                                        .collectEntries { gene_id, tx_annos ->
                                            def default_anno = tx_annos.values().first()

                                            if ( tx_annos.values().findAll { it != default_anno }.size() > 0 ) {
                                                return [ gene_id, ( tx_annos + [ 'default': 'Differing%20isoform%20descriptions' ] ) ]
                                            }

                                            [ gene_id, ( tx_annos + [ 'default': default_anno ] ) ]
                                        }

                                    def gff_lines = gff.readLines()
                                        .collect { line ->

                                            if ( line.startsWith('#') || line == '' ) { return line }

                                            def cols    = line.split('\t')
                                            def feat    = cols[2]
                                            def atts    = cols[8]

                                            if ( ! ( feat == 'gene' || feat == 'transcript' || feat == 'mRNA' ) ) { return line }

                                            def id      = feat == 'gene' ? ( atts =~ /ID=([^;]*)/ )[0][1] : ( atts =~ /Parent=([^;]*)/ )[0][1]

                                            if ( ! gene_tx_annotations.containsKey(id) ) { return line }

                                            def tx_id   = feat == 'gene' ? null : ( atts =~ /ID=([^;]*)/ )[0][1]
                                            def desc    = feat == 'gene' ? gene_tx_annotations[id]['default'] : gene_tx_annotations[id][tx_id]

                                            return ( line + ";description=$desc" )
                                        }

                                    [ "${meta.id}.described.gff" ] + gff_lines.join('\n')
                                }
                                | collectFile(newLine: true)
                                | map { file ->
                                    [ [ id: file.baseName.replace('.described', '') ], file ]
                                }

    // MODULE: GT_GFF3 as FINAL_GFF_CHECK
    ch_final_check_input        = val_describe_gff
                                ? ch_described_gff
                                : ch_target_gff

    FINAL_GFF_CHECK ( ch_final_check_input )

    ch_final_gff                = FINAL_GFF_CHECK.out.gt_gff3
    ch_versions                 = ch_versions.mix(FINAL_GFF_CHECK.out.versions.first())

    // MODULE: GFFREAD as EXTRACT_PROTEINS
    ch_extraction_inputs        = ch_final_gff
                                | join(ch_fasta)

    EXTRACT_PROTEINS(
        ch_extraction_inputs.map { meta, gff, fasta -> [ meta, gff ] },
        ch_extraction_inputs.map { meta, gff, fasta -> fasta }
    )

    ch_final_proteins           = EXTRACT_PROTEINS.out.gffread_fasta
    ch_versions                 = ch_versions.mix(EXTRACT_PROTEINS.out.versions.first())

    // MODULE: GFFREAD as EXTRACT_CDS
    ch_cds_extraction_inputs    = ch_final_gff
                                | join(ch_fasta)

    EXTRACT_CDS(
        ch_cds_extraction_inputs.map { meta, gff, fasta -> [ meta, gff ] },
        ch_cds_extraction_inputs.map { meta, gff, fasta -> fasta }
    )

    ch_final_cds                = EXTRACT_CDS.out.gffread_fasta
    ch_versions                 = ch_versions.mix(EXTRACT_CDS.out.versions.first())

    // MODULE: GFFREAD as EXTRACT_CDNA
    ch_cdna_extraction_inputs    = ch_final_gff
                                | join(ch_fasta)

    EXTRACT_CDNA(
        ch_cdna_extraction_inputs.map { meta, gff, fasta -> [ meta, gff ] },
        ch_cdna_extraction_inputs.map { meta, gff, fasta -> fasta}
    )

    ch_final_cdna                = EXTRACT_CDNA.out.gffread_fasta
    ch_versions                  = ch_versions.mix(EXTRACT_CDNA.out.versions.first())

    emit:
    final_gff                   = ch_final_gff          // [ meta, gff ]
    final_proteins              = ch_final_proteins     // [ meta, fasta ]
    final_cds                   = ch_final_cds          // [ meta, fasta ]
    final_cdna                  = ch_final_cdna         // [ meta, fasta ]
    versions                    = ch_versions           // [ versions.yml ]
}
