include { GUNZIP as GUNZIP_FASTA    } from '../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF      } from '../../modules/nf-core/gunzip'
include { GFFREAD                   } from '../../modules/nf-core/gffread'
include { LIFTOFF                   } from '../../modules/local/liftoff'

workflow FASTA_LIFTOFF {
    take:
    target_assemby              // Channel: [ meta, fasta ]
    xref_annotations_fasta      // Channel: [ meta2, fasta ]
    xref_annotations_gff        // Channel: [ meta2, gff3 ]
    
    main:
    // MODULE: GUNZIP_FASTA
    xref_annotations_fasta
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { xref_annotations_fasta_branch }

    GUNZIP_FASTA(
        xref_annotations_fasta_branch.gz
    )
    .gunzip
    | mix(
        xref_annotations_fasta_branch.rest
    )
    | set { ch_xref_annotations_gunzip_fasta }

    // MODULE: GUNZIP_GFF
    xref_annotations_gff
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { xref_annotations_gff_branch }

    GUNZIP_GFF(
        xref_annotations_gff_branch.gz
    )
    .gunzip
    | mix(
        xref_annotations_gff_branch.rest
    )
    | set { ch_xref_annotations_gunzip_gff }

    // MODULE: GFFREAD
    GFFREAD(
        ch_xref_annotations_gunzip_gff
    )
    .gff
    | set { ch_gffread_gff }

    // MODULE: LIFTOFF
    target_assemby
    | combine(
        ch_xref_annotations_gunzip_fasta
        | join(
            ch_gffread_gff
        )
    )
    | map { meta, targetFasta, refMeta, refFasta, refGFF  ->
        [[id:"${meta.id}.from.${refMeta.id}", target_assemby: meta.id], targetFasta, refFasta, refGFF]
    }
    | set { ch_liftoff_inputs }

    LIFTOFF(
        ch_liftoff_inputs.map { meta, targetFasta, refFasta, refGFF -> [meta, targetFasta] },
        ch_liftoff_inputs.map { meta, targetFasta, refFasta, refGFF -> refFasta },
        ch_liftoff_inputs.map { meta, targetFasta, refFasta, refGFF -> refGFF }
    )
    .gff3
    | map { meta, gff -> [[id: meta.target_assemby], gff] }
    | groupTuple
    | set { ch_liftoff_gff3 }

    Channel.empty()
    | mix(GUNZIP_FASTA.out.versions.first())
    | mix(GUNZIP_GFF.out.versions.first())
    | mix(GFFREAD.out.versions.first())
    | mix(LIFTOFF.out.versions.first())
    | set { ch_versions }

    emit:
    gff3        = ch_liftoff_gff3               // [ meta, [ gff3 ] ]
    versions    = ch_versions                   // [ versions.yml ]
}