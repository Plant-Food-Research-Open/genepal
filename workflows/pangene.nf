include { validateParams                        } from '../modules/local/validate_params'
include { id_from_file_name                     } from '../modules/local/validate_params'
include { PREPARE_ASSEMBLY                      } from '../subworkflows/local/prepare_assembly'
include { PREPROCESS_RNASEQ                     } from '../subworkflows/local/preprocess_rnaseq'
include { ALIGN_RNASEQ                          } from '../subworkflows/local/align_rnaseq'
include { PREPARE_EXT_PROTS                     } from '../subworkflows/local/prepare_ext_prots'
include { FASTA_BRAKER3                         } from '../subworkflows/local/fasta_braker3'
include { FASTA_LIFTOFF                         } from '../subworkflows/local/fasta_liftoff'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions'

validateParams(params)

workflow PANGENE {

    ch_versions                 = Channel.empty()

    ch_target_assembly          = Channel.fromList(params.target_assemblies)
                                | map { it ->
                                    def tag     = it[0]
                                    def fasta   = it[1]

                                    [ [ id: tag ], file(fasta, checkIfExists: true) ]
                                }

    ch_braker_annotation        = Channel.fromList(params.target_assemblies)
                                | map { it ->
                                    if ( it.size() == 4 ) {
                                        def tag         = it[0]
                                        def braker_gff3 = it[2]
                                        def hints_gff   = it[3]

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

    ch_te_library               = Channel.fromList(params.te_libraries)
                                | map { tag, fasta ->
                                    [ [ id:tag ], file(fasta, checkIfExists: true) ]
                                }

    ch_samplesheet              = params.samplesheet
                                ? Channel.fromPath(params.samplesheet, checkIfExists: true)
                                : Channel.empty()

    ch_tar_assm_str             = Channel.of(
                                    params.target_assemblies
                                    .collect { it -> it[0].strip() }.join(",") // it[0] = tag
                                )

    ch_ribo_db                  = params.remove_ribo_rna
                                ? file(params.ribo_database_manifest, checkIfExists: true)
                                : null

    ch_sortmerna_fastas         = ch_ribo_db
                                ? Channel.from(ch_ribo_db ? ch_ribo_db.readLines() : null)
                                | map { row -> file(row, checkIfExists: true) }
                                | collect
                                : Channel.empty()

    ch_ext_prot_fastas          = params.external_protein_fastas
                                ? Channel.fromList(params.external_protein_fastas)
                                | map { filePath ->
                                    def fileHandle = file(filePath, checkIfExists: true)
                                    [ [ id: id_from_file_name( fileHandle.baseName ) ], fileHandle]
                                }
                                : Channel.empty()

    ch_xref_mm                  = params.liftoff_xref_annotations
                                ? Channel.fromList(params.liftoff_xref_annotations)
                                | multiMap { fasta, gff ->
                                    def fastaFile = file(fasta, checkIfExists:true)

                                    fasta: [ [ id: id_from_file_name( fastaFile.baseName ) ], fastaFile ]
                                    gff: [ [ id: id_from_file_name( fastaFile.baseName ) ], file(gff, checkIfExists:true) ]
                                }
                                : Channel.empty()

    ch_xref_fasta               = ch_xref_mm.fasta
    ch_xref_gff                 = ch_xref_mm.gff

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
        ch_samplesheet,
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
        ch_xref_fasta,
        ch_xref_gff
    )

    ch_liftoff_gff3             = FASTA_LIFTOFF.out.gff3
    ch_versions                 = ch_versions.mix(FASTA_LIFTOFF.out.versions)

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}
