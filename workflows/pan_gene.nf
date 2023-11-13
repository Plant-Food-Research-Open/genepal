nextflow.enable.dsl=2

include { BRAKER3                               } from '../modules/kherronism/braker3'
include { GUNZIP as GUNZIP_XREF_FASTA           } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_XREF_GFF             } from '../modules/nf-core/gunzip'
include { validateParams                        } from '../modules/local/validate_params'

include { PREPARE_ASSEMBLY                      } from '../subworkflows/local/prepare_assembly'
include { PREPROCESS_RNASEQ                     } from '../subworkflows/local/preprocess_rnaseq'
include { ALIGN_RNASEQ                          } from '../subworkflows/local/align_rnaseq'
include { PREPARE_EXT_PROTS                     } from '../subworkflows/local/prepare_ext_prots'

validateParams(params)

workflow PAN_GENE {

    ch_versions                 = Channel.empty()

    ch_target_assembly          = Channel.fromList(params.target_assemblies)
                                | map { tag, filePath ->
                                    [[id:tag], file(filePath, checkIfExists: true)]
                                }

    ch_te_library               = Channel.fromList(params.te_libraries)
                                | map { tag, filePath ->
                                    [[id:tag], file(filePath, checkIfExists: true)]
                                }

    ch_samplesheet              = params.samplesheet
                                ? Channel.fromPath(params.samplesheet, checkIfExists: true)
                                : Channel.empty()
    
    ch_tar_assm_str             = Channel.of(
                                    params.target_assemblies
                                    .collect { tag, fastaPath -> tag.strip() }.join(",")
                                )

    ch_ribo_db                  = params.remove_ribo_rna
                                ? file(params.ribo_database_manifest, checkIfExists: true)
                                : null

    ch_sortmerna_fastas         = Channel.from(ch_ribo_db ? ch_ribo_db.readLines() : null)
                                | map { row -> file(row, checkIfExists: true) }
                                | collect

    ch_ext_prot_fastas          = (params.external_protein_fastas
                                ? Channel.fromList(params.external_protein_fastas)
                                : Channel.empty())
                                | map { filePath ->
                                    def fileHandle = file(filePath, checkIfExists: true)
                                    [[id:fileHandle.getSimpleName()], fileHandle]
                                }

    // SUBWORKFLOW: PREPARE_ASSEMBLY
    PREPARE_ASSEMBLY(
        ch_target_assembly,
        ch_te_library
    )

    ch_valid_target_assembly    = PREPARE_ASSEMBLY.out.target_assemby
    ch_masked_target_assembly   = PREPARE_ASSEMBLY.out.masked_target_assembly
    ch_target_assemby_index     = PREPARE_ASSEMBLY.out.target_assemby_index
    ch_versions                 = ch_versions.mix(PREPARE_ASSEMBLY.out.versions)

    // SUBWORKFLOW: PREPROCESS_RNASEQ
    PREPROCESS_RNASEQ(
        ch_samplesheet,
        ch_tar_assm_str,
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
        ch_target_assemby_index
    )

    ch_rnaseq_bam               = ALIGN_RNASEQ.out.bam
    ch_versions                 = ch_versions.mix(ALIGN_RNASEQ.out.versions)

    // MODULE: PREPARE_EXT_PROTS
    PREPARE_EXT_PROTS(
        ch_ext_prot_fastas
    )

    ch_ext_prots_fasta          = PREPARE_EXT_PROTS.out.ext_prots_fasta
    ch_versions                 = ch_versions.mix(PREPARE_EXT_PROTS.out.versions)

    // MODULE: BRAKER3
    ch_braker_inputs            = ch_masked_target_assembly
                                | join(ch_rnaseq_bam, remainder: true)
                                | combine(
                                    ch_ext_prots_fasta.map { meta, filePath -> filePath }.ifEmpty(null)
                                )
                                | map { meta, fasta, bam, prots -> [meta, fasta, bam ?: [], prots ?: []] }
    
    def rnaseq_sets_dirs        = []
    def rnaseq_sets_ids         = []
    def hintsfile               = []

    BRAKER3(
        ch_braker_inputs.map { meta, fasta, bam, prots -> [meta, fasta] },
        ch_braker_inputs.map { meta, fasta, bam, prots -> bam },
        rnaseq_sets_dirs,
        rnaseq_sets_ids,
        ch_braker_inputs.map { meta, fasta, bam, prots -> prots },
        hintsfile
    )

    ch_braker_gff3              = BRAKER3.out.gff3
    ch_versions                 = ch_versions.mix(BRAKER3.out.versions.first())

    // // MODULE: GUNZIP_XREF_FASTA
    // ch_xref_annotations = Channel.empty()
    // if(params.liftoff_xref_annotations) {
    //     Channel.fromList(params.liftoff_xref_annotations)
    //     | multiMap { fasta, gff ->
    //         def fastaFile = file(fasta, checkIfExists:true)
    //         def meta = [id:fastaFile.getSimpleName()]

    //         fasta: [meta, fastaFile]
    //         gff: [meta, file(gff, checkIfExists:true)]
    //     }
    //     | set { ch_xref_annotations }
    // }

    // ch_xref_annotations.fasta
    // | branch { meta, file ->
    //     gz: "$file".endsWith(".gz")
    //     rest: !"$file".endsWith(".gz")
    // }
    // | set { ch_xref_annotations_branch }

    // GUNZIP_XREF_FASTA(
    //     ch_xref_annotations_branch.gz
    // )
    // .gunzip
    // | mix(
    //     ch_xref_annotations_branch.rest
    // )
    // | set { ch_xref_annotations_fasta }

    // // MODULE: GUNZIP_XREF_GFF
    // ch_xref_annotations.gff
    // | branch { meta, file ->
    //     gz: "$file".endsWith(".gz")
    //     rest: !"$file".endsWith(".gz")
    // }
    // | set { ch_xref_annotations_gff_branch }

    // GUNZIP_XREF_GFF(
    //     ch_xref_annotations_gff_branch.gff.map { meta, fasta, gff -> [meta, gff] }
    // )
    // .gunzip
    // | mix(
    //     ch_xref_annotations_gff_branch.rest.map { meta, fasta, gff -> [meta, gff] }
    // )
    // | set { ch_xref_annotations_gff }

    // ch_xref_annotations_fasta
    // | join(
    //     ch_xref_annotations_gff
    // )
    // | set { ch_xref_annotations }

    // // MODULE: LIFTOFF
    // ch_xref_annotations
    // | combine(
    //     ch_validated_target_assemblies
    // )
    // | map { meta, ref_fasta, refGFF, targetMeta, targetFasta -> [[id:"${targetMeta.id}.from.${meta.id}"], ref_fasta, refGFF, targetFasta] }
    // | set { ch_liftoff_inputs }
}