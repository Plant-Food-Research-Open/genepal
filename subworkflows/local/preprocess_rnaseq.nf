include { CAT_FASTQ                     } from '../../modules/nf-core/cat/fastq'
include { SORTMERNA                     } from '../../modules/nf-core/sortmerna'
include { EXTRACT_SAMPLES               } from '../../subworkflows/local/extract_samples'
include { FASTQ_FASTQC_UMITOOLS_FASTP   } from '../../subworkflows/nf-core/fastq_fastqc_umitools_fastp'

workflow PREPROCESS_RNASEQ {
    take:
    samplesheet                     // path: csv
    permissible_target_assemblies   // val: assembly_a,assembly_b
    skip_fastqc                     // val: true|false
    skip_fastp                      // val: true|false
    save_trimmed                    // val: true|false
    min_trimmed_reads               // val: Integer
    remove_ribo_rna                 // val: true|false
    sortmerna_fastas                // channel: [ [ fasta ] ]
    
    main:
    ch_versions = Channel.empty()
    // SUBWORKFLOW: EXTRACT_SAMPLES
    EXTRACT_SAMPLES(
        samplesheet,
        ch_permissible_target_assemblies
    )
    .reads
    | map { meta, fastq ->
        groupID = meta.id - ~/_T\d+/
        [ meta + [id: groupID], fastq ]
    }
    | groupTuple()
    | branch { meta, fastq ->
        single  : fastq.size() == 1
            return [ meta, fastq.flatten() ]
        multiple: fastq.size() > 1
            return [ meta, fastq.flatten() ]
    }
    | set { ch_fastq }

    EXTRACT_SAMPLES.out.assemblies
    | map { meta, assembly ->
        groupID = meta.id - ~/_T\d+/
        [ meta + [id: groupID], assembly ]
    }
    | unique
    | set { ch_reads_target }

    // MODULES: CAT_FASTQ
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    | mix(ch_fastq.single)
    | set { ch_cat_fastq }

    // SUBWORKFLOW: FASTQ_FASTQC_UMITOOLS_FASTP
    def with_umi            = false
    def skip_umi_extract    = true
    def umi_discard_read    = false
    FASTQ_FASTQC_UMITOOLS_FASTP (
        ch_cat_fastq,
        skip_fastqc,
        with_umi,
        skip_umi_extract,
        umi_discard_read,
        skip_fastp,
        [],
        save_trimmed,
        save_trimmed,
        min_trimmed_reads
    )
    .reads
    | set { ch_trim_reads }

    // MODULE: SORTMERNA
    if (remove_ribo_rna) {
        SORTMERNA (
            ch_trim_reads,
            sortmerna_fastas
        )
        .reads
        | set { ch_sortmerna_reads }

        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
    }

    ch_versions
    | mix(EXTRACT_SAMPLES.out.versions)
    | mix(CAT_FASTQ.out.versions.first())
    | mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)
    | set { ch_versions }

    emit:
    trim_reads      = remove_ribo_rna ? ch_sortmerna_reads : ch_trim_reads  // channel: [ meta, [ fq ] ]
    reads_target    = ch_reads_target                                       // channel: [ meta, assembly_id ]
    versions        = ch_versions                                           // channel: [ versions.yml ]
}