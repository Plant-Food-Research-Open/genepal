nextflow.enable.dsl=2


include { GUNZIP as GUNZIP_EXTERNAL_PROTEIN_SEQ } from '../modules/nf-core/gunzip'
include { CAT_CAT as CAT_PROTEIN_SEQS           } from '../modules/nf-core/cat/cat'
include { BRAKER3                               } from '../modules/kherronism/braker3'
include { GUNZIP as GUNZIP_XREF_FASTA           } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_XREF_GFF             } from '../modules/nf-core/gunzip'
include { validateParams                        } from '../modules/local/validate_params'

include { PREPARE_ASSEMBLY                      } from '../subworkflows/local/prepare_assembly'
include { PREPROCESS_RNASEQ                     } from '../subworkflows/local/preprocess_rnaseq'
include { ALIGN_RNASEQ                          } from '../subworkflows/local/align_rnaseq'

validateParams(params)

// Additional validation
// Check rRNA databases for sortmerna
if (params.sample_prep.remove_ribo_rna) {
    ch_ribo_db = file(params.sample_prep.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

workflow PAN_GENE {

    // Versions
    ch_versions = Channel.empty()

    // Input channels
    Channel.fromList(params.target_assemblies)
    | map { tag, filePath ->
        [[id:tag], file(filePath, checkIfExists: true)]
    }
    | set { ch_target_assembly }

    Channel.fromList(params.te_libraries)
    | map { tag, filePath ->
        [[id:tag], file(filePath, checkIfExists: true)]
    }
    | set { ch_te_library }

    ch_samplesheet = Channel.empty()
    if(params.samplesheet) {
        ch_samplesheet = Channel.fromPath(params.samplesheet)
    }
    
    Channel.of(params.target_assemblies.collect { tag, fastaPath -> tag.strip() }.join(","))
    | set { ch_permissible_target_assemblies }

    Channel.from(ch_ribo_db.readLines())
    | map { row -> file(row, checkIfExists: true) }
    | collect
    | set { ch_sortmerna_fastas }

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
        ch_permissible_target_assemblies,
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

    // SUBWORKFLOW: STAR_ALIGN
    ALIGN_RNASEQ(
        ch_reads_target,
        ch_trim_reads,
        ch_target_assemby_index
    )

    // MODULE: GUNZIP_EXTERNAL_PROTEIN_SEQ
    ch_ext_prot_seqs = Channel.empty()
    if(params.external_protein_seqs) {
        ch_ext_prot_seqs = Channel.fromList(params.external_protein_seqs)
    }
    
    ch_ext_prot_seqs
    | map { filePath ->
        def fileHandle = file(filePath, checkIfExists: true)
        [[id:fileHandle.getSimpleName()], fileHandle]
    }
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_ext_prot_seqs_branch }

    GUNZIP_EXTERNAL_PROTEIN_SEQ(
        ch_ext_prot_seqs_branch.gz
    )
    .gunzip
    | mix(
        ch_ext_prot_seqs_branch.rest
    )
    | set { ch_ext_prot_seqs }

    ch_versions = ch_versions.mix(GUNZIP_EXTERNAL_PROTEIN_SEQ.out.versions.first())

    // MODULE: CAT_PROTEIN_SEQS
    ch_ext_prot_seqs
    | map { meta, filePath -> filePath }
    | collect
    | map { fileList -> [[id:"protein_seqs"], fileList] }
    | CAT_PROTEIN_SEQS
    
    ch_ext_prot_seqs = CAT_PROTEIN_SEQS.out.file_out
    ch_versions = ch_versions.mix(CAT_PROTEIN_SEQS.out.versions)

    // MODULE: BRAKER3
    REPEATMASKER.out.fasta_masked
    | mix(ch_samtools_bam)
    | groupTuple(size: 2, remainder: true)
    | map { meta, groupedItems ->
        def maskedFasta = groupedItems[0]

        if(groupedItems.size() == 2) {
            def bam = groupedItems[1]
            return [meta, maskedFasta, bam]
        } else {
            return [meta, maskedFasta, []]
        }
    }
    | set { ch_braker_inputs }
    
    if(params.external_protein_seqs) {
        ch_braker_inputs
        | combine(ch_ext_prot_seqs.map{meta, filePath -> filePath})
        | set { ch_braker_inputs }
    } else {
        ch_braker_inputs
        | map{meta, assembly, bam -> [meta, assembly, bam, []]}
        | set { ch_braker_inputs }
    }
    
    ch_fasta            = ch_braker_inputs.map { meta, assembly, bam, proteinSeq -> [meta, assembly] }
    ch_bam              = ch_braker_inputs.map { meta, assembly, bam, proteinSeq -> bam }
    ch_proteins         = ch_braker_inputs.map { meta, assembly, bam, proteinSeq -> proteinSeq }
    ch_rnaseq_sets_dirs = []
    ch_rnaseq_sets_ids  = []
    ch_hintsfile        = []

    BRAKER3(
        ch_fasta,
        ch_bam,
        ch_rnaseq_sets_dirs,
        ch_rnaseq_sets_ids,
        ch_proteins,
        ch_hintsfile
    )

    ch_versions = ch_versions.mix(BRAKER3.out.versions.first())

    // MODULE: GUNZIP_XREF_FASTA
    ch_xref_annotations = Channel.empty()
    if(params.liftoff.xref_annotations) {
        Channel.fromList(params.liftoff.xref_annotations)
        | multiMap { fasta, gff ->
            def fastaFile = file(fasta, checkIfExists:true)
            def meta = [id:fastaFile.getSimpleName()]

            fasta: [meta, fastaFile]
            gff: [meta, file(gff, checkIfExists:true)]
        }
        | set { ch_xref_annotations }
    }

    ch_xref_annotations.fasta
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_xref_annotations_branch }

    GUNZIP_XREF_FASTA(
        ch_xref_annotations_branch.gz
    )
    .gunzip
    | mix(
        ch_xref_annotations_branch.rest
    )
    | set { ch_xref_annotations_fasta }

    // MODULE: GUNZIP_XREF_GFF
    ch_xref_annotations.gff
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_xref_annotations_gff_branch }

    GUNZIP_XREF_GFF(
        ch_xref_annotations_gff_branch.gff.map { meta, fasta, gff -> [meta, gff] }
    )
    .gunzip
    | mix(
        ch_xref_annotations_gff_branch.rest.map { meta, fasta, gff -> [meta, gff] }
    )
    | set { ch_xref_annotations_gff }

    ch_xref_annotations_fasta
    | join(
        ch_xref_annotations_gff
    )
    | set { ch_xref_annotations }

    // // MODULE: LIFTOFF
    // ch_xref_annotations
    // | combine(
    //     ch_validated_target_assemblies
    // )
    // | map { meta, ref_fasta, refGFF, targetMeta, targetFasta -> [[id:"${targetMeta.id}.from.${meta.id}"], ref_fasta, refGFF, targetFasta] }
    // | set { ch_liftoff_inputs }
}