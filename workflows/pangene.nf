include { fromSamplesheet; paramsSummaryLog     } from 'plugin/nf-validation'
include { idFromFileName; validateFastqMetadata } from '../modules/local/utils'
include { PREPARE_ASSEMBLY                      } from '../subworkflows/local/prepare_assembly'
include { PREPROCESS_RNASEQ                     } from '../subworkflows/local/preprocess_rnaseq'
include { ALIGN_RNASEQ                          } from '../subworkflows/local/align_rnaseq'
include { PREPARE_EXT_PROTS                     } from '../subworkflows/local/prepare_ext_prots'
include { FASTA_BRAKER3                         } from '../subworkflows/local/fasta_braker3'
include { FASTA_LIFTOFF                         } from '../subworkflows/local/fasta_liftoff'
include { MERGE_ANNOTATIONS                     } from '../subworkflows/local/merge_annotations'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions'

log.info paramsSummaryLog(workflow)

workflow PANGENE {

    // Versions channel
    ch_versions                 = Channel.empty()

    // Input channels
    ch_input                    = Channel.fromSamplesheet('input')

    ch_target_assembly          = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def fasta       = it[1]

                                    [ [ id: tag ], file(fasta, checkIfExists: true) ]
                                }

    ch_tar_assm_str             = ch_input
                                | map { it ->
                                    def tag         = it[0].strip()

                                    tag
                                }
                                | collect
                                | map { it ->
                                    it.join(",")
                                }

    ch_masked                   = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def is_masked   = it[2]

                                    [ [ id: tag ], is_masked == "yes" ]
                                }

    ch_te_library               = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def te_fasta    = it[3]

                                    if ( te_fasta ) {
                                        [ [ id:tag ], file(te_fasta, checkIfExists: true) ]
                                    }
                                }

    ch_braker_annotation        = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def braker_gff3 = it[4]
                                    def hints_gff   = it[5]

                                    if ( braker_gff3 ) {
                                        [
                                            [ id: tag ],
                                            file(braker_gff3, checkIfExists: true),
                                            file(hints_gff, checkIfExists: true)
                                        ]
                                    }
                                }

    ch_braker_ex_asm_str        = ch_braker_annotation
                                | map { meta, braker_gff3, hints_gff -> meta.id }
                                | collect
                                | map { it.join(",") }
                                | ifEmpty( "" )

    ch_reads                    = ! params.fastq
                                ? Channel.empty()
                                : Channel.fromSamplesheet('fastq')
                                | map { meta, fq1, fq2 ->
                                    fq2
                                    ? [ meta + [ single_end: false ], [ file(fq1, checkIfExists:true), file(fq2, checkIfExists:true) ] ]
                                    : [ meta + [ single_end: true ], [ file(fq1, checkIfExists:true) ] ]
                                }
                                | map { meta, fqs ->
                                    [ meta.id, meta + [ target_assemblies: meta.target_assemblies.split(';').sort() ], fqs ]
                                }
                                | groupTuple
                                | combine(ch_tar_assm_str)
                                | map { id, metas, fqs, tar_assm_str ->
                                    validateFastqMetadata(metas, fqs, tar_assm_str)
                                }

    ch_ribo_db                  = params.remove_ribo_rna
                                ? file(params.ribo_database_manifest, checkIfExists: true)
                                : null

    ch_sortmerna_fastas         = ch_ribo_db
                                ? Channel.from(ch_ribo_db ? ch_ribo_db.readLines() : null)
                                | map { row -> file(row, checkIfExists: true) }
                                | collect
                                : Channel.empty()

    ch_ext_prot_fastas          = ! params.external_protein_fastas
                                ? Channel.empty()
                                : Channel.fromPath(params.external_protein_fastas)
                                | splitText
                                | map { file_path ->
                                    def file_handle = file(file_path.strip(), checkIfExists: true)
                                    [ [ id: idFromFileName( file_handle.baseName ) ], file_handle ]
                                }

    ch_liftoff_mm               = ! params.liftoff_annotations
                                ? Channel.empty()
                                : Channel.fromSamplesheet('liftoff_annotations')
                                | multiMap { fasta, gff ->
                                    def fastaFile = file(fasta, checkIfExists:true)

                                    fasta: [ [ id: idFromFileName( fastaFile.baseName ) ], fastaFile ]
                                    gff: [ [ id: idFromFileName( fastaFile.baseName ) ], file(gff, checkIfExists:true) ]
                                }

    ch_liftoff_fasta               = ch_liftoff_mm.fasta
    ch_liftoff_gff                 = ch_liftoff_mm.gff

    // SUBWORKFLOW: PREPARE_ASSEMBLY
    PREPARE_ASSEMBLY(
        ch_target_assembly,
        ch_te_library,
        params.repeat_annotator,
        ch_braker_ex_asm_str
    )

    ch_valid_target_assembly    = PREPARE_ASSEMBLY.out.target_assemby
    ch_masked_target_assembly   = PREPARE_ASSEMBLY.out.masked_target_assembly
    ch_target_assemby_index     = PREPARE_ASSEMBLY.out.target_assemby_index
    ch_versions                 = ch_versions.mix(PREPARE_ASSEMBLY.out.versions)

    // SUBWORKFLOW: PREPROCESS_RNASEQ
    PREPROCESS_RNASEQ(
        ch_reads,
        ch_tar_assm_str,
        ch_braker_ex_asm_str,
        params.skip_fastqc,
        params.skip_fastp,
        params.save_trimmed,
        params.min_trimmed_reads,
        params.remove_ribo_rna,
        ch_sortmerna_fastas
    )

    ch_trim_reads               = PREPROCESS_RNASEQ.out.trim_reads
    ch_reads_target             = PREPROCESS_RNASEQ.out.reads_target
    ch_versions                 = ch_versions.mix(PREPROCESS_RNASEQ.out.versions)

    // SUBWORKFLOW: ALIGN_RNASEQ
    ALIGN_RNASEQ(
        ch_reads_target,
        ch_trim_reads,
        ch_target_assemby_index,
    )

    ch_rnaseq_bam               = ALIGN_RNASEQ.out.bam
    ch_versions                 = ch_versions.mix(ALIGN_RNASEQ.out.versions)

    // MODULE: PREPARE_EXT_PROTS
    PREPARE_EXT_PROTS(
        ch_ext_prot_fastas
    )

    ch_ext_prots_fasta          = PREPARE_EXT_PROTS.out.ext_prots_fasta
    ch_versions                 = ch_versions.mix(PREPARE_EXT_PROTS.out.versions)

    // SUBWORKFLOW: FASTA_BRAKER3
    FASTA_BRAKER3(
        ch_masked_target_assembly,
        ch_braker_ex_asm_str,
        ch_rnaseq_bam,
        ch_ext_prots_fasta,
        ch_braker_annotation
    )

    ch_braker_gff3              = FASTA_BRAKER3.out.braker_gff3
    ch_braker_hints             = FASTA_BRAKER3.out.braker_hints
    ch_versions                 = ch_versions.mix(FASTA_BRAKER3.out.versions)

    // SUBWORKFLOW: FASTA_LIFTOFF
    FASTA_LIFTOFF(
        ch_valid_target_assembly,
        ch_liftoff_fasta,
        ch_liftoff_gff
    )

    ch_liftoff_gff3             = FASTA_LIFTOFF.out.gff3
    ch_versions                 = ch_versions.mix(FASTA_LIFTOFF.out.versions)

    // SUBWORKFLOW: MERGE_ANNOTATIONS
    MERGE_ANNOTATIONS ( ch_braker_gff3, ch_liftoff_gff3 )

    ch_merged_gff               = MERGE_ANNOTATIONS.out.merged_gff
    ch_versions                 = ch_versions.mix(MERGE_ANNOTATIONS.out.versions)

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
