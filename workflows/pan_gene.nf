nextflow.enable.dsl=2

include { GUNZIP as GUNZIP_TARGET_ASSEMBLY      } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_TE_LIBRARY           } from '../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_EXTERNAL_PROTEIN_SEQ } from '../modules/nf-core/gunzip'
include { FASTA_VALIDATE                        } from '../modules/local/fasta_validate'
include { REPEATMASKER                          } from '../modules/kherronism/repeatmasker'
include { STAR_GENOMEGENERATE                   } from '../modules/nf-core/star/genomegenerate'
include { CAT_FASTQ                             } from '../modules/nf-core/cat/fastq'
include { SORTMERNA                             } from '../modules/nf-core/sortmerna'
include { STAR_ALIGN                            } from '../modules/nf-core/star/align'
include { SAMTOOLS_CAT                          } from '../modules/nf-core/samtools/cat'
include { CAT_CAT as CAT_PROTEIN_SEQS           } from '../modules/nf-core/cat/cat'
include { BRAKER3                               } from '../modules/kherronism/braker3'

include { PERFORM_EDTA_ANNOTATION               } from '../subworkflows/local/perform_edta_annotation'
include { EXTRACT_SAMPLES                       } from '../subworkflows/local/extract_samples'
include { FASTQ_FASTQC_UMITOOLS_FASTP           } from '../subworkflows/nf-core/fastq_fastqc_umitools_fastp'

include { validateParams                        } from '../modules/local/validate_params'

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
    
    // MODULE: GUNZIP_TARGET_ASSEMBLY
    Channel.fromList(params.target_assemblies)
    | map { tag, filePath ->
        [[id:tag], file(filePath, checkIfExists: true)]
    }
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_target_assemblies }

    GUNZIP_TARGET_ASSEMBLY(
        ch_target_assemblies.gz
    )
    .gunzip
    | mix(
        ch_target_assemblies.rest
    )
    | set { ch_gunzip_target_assemblies }

    ch_versions = ch_versions.mix(GUNZIP_TARGET_ASSEMBLY.out.versions.first())

    // MODULE: FASTA_VALIDATE
    FASTA_VALIDATE(ch_gunzip_target_assemblies)
    .valid_fasta
    | set { ch_validated_target_assemblies }

    ch_versions = ch_versions.mix(FASTA_VALIDATE.out.versions.first())

    // MODULE: GUNZIP_TE_LIBRARY
    Channel.fromList(params.te_libraries)
    | map { tag, filePath ->
        [[id:tag], file(filePath, checkIfExists: true)]
    }
    | branch { meta, file ->
        gz: "$file".endsWith(".gz")
        rest: !"$file".endsWith(".gz")
    }
    | set { ch_te_libraries }

    GUNZIP_TE_LIBRARY(
        ch_te_libraries.gz
    )
    .gunzip
    | mix(
        ch_te_libraries.rest
    )
    | set { ch_gunzip_te_libraries }

    ch_versions = ch_versions.mix(GUNZIP_TE_LIBRARY.out.versions.first())

    // SUBWORKFLOW: PERFORM_EDTA_ANNOTATION
    ch_validated_target_assemblies
    | join(
        ch_gunzip_te_libraries, remainder: true
    )
    | filter { meta, assembly, teLib ->
        teLib == null
    }
    | map {meta, assembly, teLib -> [meta, assembly]}
    | PERFORM_EDTA_ANNOTATION

    ch_versions = ch_versions.mix(PERFORM_EDTA_ANNOTATION.out.versions)
    
    // MODULE: REPEATMASKER
    ch_validated_target_assemblies
    | join(
        PERFORM_EDTA_ANNOTATION.out.te_lib_fasta.mix(ch_gunzip_te_libraries)
    )
    | set { ch_assemblies_n_te_libs }

    REPEATMASKER(
        ch_assemblies_n_te_libs.map {meta, assembly, teLib -> [meta, assembly]},
        ch_assemblies_n_te_libs.map {meta, assembly, teLib -> teLib},
    )

    ch_versions = ch_versions.mix(REPEATMASKER.out.versions.first())

    // MODULE: STAR_GENOMEGENERATE
    def star_ignore_sjdbgtf = true
    STAR_GENOMEGENERATE(
        REPEATMASKER.out.fasta_masked,
        REPEATMASKER.out.fasta_masked.map{meta, maskedFasta -> [meta, []]},
        star_ignore_sjdbgtf
    )
    .index
    | set { ch_assembly_index }

    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

    // SUBWORKFLOW: EXTRACT_SAMPLES
    ch_samplesheet_path = Channel.empty()
    if(params.samplesheet != null) {
        ch_samplesheet_path = Channel.fromPath(params.samplesheet)
    }
    
    EXTRACT_SAMPLES(
        ch_samplesheet_path,
        Channel.of(params.target_assemblies.collect{tag, fastaPath -> tag.strip()}.join(","))
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

    ch_read_target_assemblies = EXTRACT_SAMPLES.out.assemblies
    ch_versions = ch_versions.mix(EXTRACT_SAMPLES.out.versions)

    // MODULES: CAT_FASTQ
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    | mix(ch_fastq.single)
    | set { ch_cat_fastq }
    
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    // SUBWORKFLOW: FASTQ_FASTQC_UMITOOLS_FASTP
    def with_umi            = false
    def skip_umi_extract    = true
    def umi_discard_read    = false
    FASTQ_FASTQC_UMITOOLS_FASTP (
        ch_cat_fastq,
        params.sample_prep.skip_fastqc,
        with_umi,
        skip_umi_extract,
        umi_discard_read,
        params.sample_prep.skip_fastp,
        [],
        params.sample_prep.save_trimmed,
        params.sample_prep.save_trimmed,
        params.sample_prep.min_trimmed_reads
    )
    .reads
    | set { ch_trim_reads }

    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)

    // MODULE: SORTMERNA
    if (params.sample_prep.remove_ribo_rna) {
        Channel.from(ch_ribo_db.readLines())
        | map { row -> file(row, checkIfExists: true) }
        | collect
        | set { ch_sortmerna_fastas }

        SORTMERNA (
            ch_trim_reads,
            ch_sortmerna_fastas
        )
        .reads
        | set { ch_trim_reads }

        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
    }

    // MODULE: STAR_ALIGN
    ch_read_target_assemblies
    | map { meta, assembly ->
        groupID = meta.id - ~/_T\d+/
        [ meta + [id: groupID], assembly ]
    }
    | unique
    | combine(ch_trim_reads, by:0)
    | map { meta, assembly, fastq ->
        [assembly, [id:"${meta.id}.on.${assembly}", single_end:meta.single_end, target_assembly:assembly], fastq]
    }
    | combine(
        ch_assembly_index.map { meta, index -> [meta.id, index] },
        by:0
    )
    | map { assembly, meta, fastq, index -> [meta, fastq, index] }
    | set { ch_star_inputs }

    def seq_platform = false
    def seq_center = false
    STAR_ALIGN(
        ch_star_inputs.map{meta, fastq, index -> [meta, fastq]},
        ch_star_inputs.map{meta, fastq, index -> [[id: meta.target_assembly], index]},
        ch_star_inputs.map{meta, fastq, index -> [[id: meta.target_assembly], []]},
        star_ignore_sjdbgtf,
        seq_platform,
        seq_center
    )
    .bam_sorted
    | set { ch_star_bam }

    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    // MODULE: SAMTOOLS_CAT
    ch_star_bam
    | map { meta, bam ->
        [
            [id: meta.target_assembly],
            bam instanceof List ? bam.find {it =~ /Aligned/} : bam
        ]
    }
    | groupTuple
    | branch { meta, bamList ->
        bams: bamList.size() > 1
        bam: bamList.size() <= 1
    }
    | set { ch_star_bam_branch }

    SAMTOOLS_CAT(
        ch_star_bam_branch.bams
    )
    .bam.map { meta, bam -> [meta, [bam]] }
    | mix(
        ch_star_bam_branch.bam
    )
    | set { ch_samtools_bam }

    ch_versions = ch_versions.mix(SAMTOOLS_CAT.out.versions.first())

    // MODULE: GUNZIP_EXTERNAL_PROTEIN_SEQ
    ch_ext_prot_seqs = Channel.empty()
    if(params.external_protein_seqs != null) {
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
    | map{ meta, filePath -> filePath }
    | collect
    | map{ fileList -> [[id:"protein_seqs"], fileList] }
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
    
    ch_fasta            = ch_braker_inputs.map{ meta, assembly, bam, proteinSeq -> [meta, assembly] }
    ch_bam              = ch_braker_inputs.map{ meta, assembly, bam, proteinSeq -> bam }
    ch_proteins         = ch_braker_inputs.map{ meta, assembly, bam, proteinSeq -> proteinSeq }
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
}