include { CAT_FASTQ                     } from '../../modules/nf-core/cat/fastq'
include { SORTMERNA                     } from '../../modules/nf-core/sortmerna'
include { EXTRACT_SAMPLES               } from '../../subworkflows/local/extract_samples'
include { FASTQ_FASTQC_UMITOOLS_FASTP   } from '../../subworkflows/nf-core/fastq_fastqc_umitools_fastp'

workflow PREPROCESS_RNASEQ {
    take:
    samplesheet                     // path: csv
    permissible_assemblies          // val: assembly_a,assembly_b
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
        permissible_assemblies
    )

    ch_fastq                        = EXTRACT_SAMPLES.out.reads
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

    ch_reads_target                 = EXTRACT_SAMPLES.out.assemblies
                                    | map { meta, assembly ->
                                        groupID = meta.id - ~/_T\d+/
                                        [ meta + [id: groupID], assembly ]
                                    }
                                    | unique

    ch_versions                     = ch_versions.mix(EXTRACT_SAMPLES.out.versions)

    // MODULES: CAT_FASTQ
    CAT_FASTQ ( ch_fastq.multiple )

    ch_cat_fastq                    = CAT_FASTQ.out.reads.mix(ch_fastq.single)
    ch_versions                     = ch_versions.mix(CAT_FASTQ.out.versions.first())

    // SUBWORKFLOW: FASTQ_FASTQC_UMITOOLS_FASTP
    def with_umi                    = false
    def skip_umi_extract            = true
    def umi_discard_read            = false

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

    ch_trim_reads                   = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads

    ch_cat_fastq
    | join(ch_trim_reads, remainder:true)
    | map { meta, reads, trimmed ->
        if (!trimmed) {
            System.err.println("WARNING: Dropping ${reads.collect { it.getName() }} as read count after trimming is less than $min_trimmed_reads")
        }
    }

    ch_versions                     = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions.first())

    // MODULE: SORTMERNA
    SORTMERNA(
        remove_ribo_rna ? ch_trim_reads : Channel.empty(),
        sortmerna_fastas.map { fastas -> [ [], fastas ] },
        [ [], [] ]
    )

    ch_emitted_reads                = remove_ribo_rna
                                    ? SORTMERNA.out.reads
                                    : ch_trim_reads
    ch_versions                     = ch_versions.mix(SORTMERNA.out.versions.first())



    emit:
    trim_reads                      = ch_emitted_reads  // channel: [ meta, [ fq ] ]
    reads_target                    = ch_reads_target   // channel: [ meta, assembly_id ]
    versions                        = ch_versions       // channel: [ versions.yml ]
}
