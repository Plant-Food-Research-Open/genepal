//
// Subworkflow with functionality specific to the plant-food-research-open/genepal pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version                 // boolean: Display version and exit
    monochrome_logs         // boolean: Do not use coloured log outputs
    nextflow_cli_args       //   array: List of positional nextflow CLI args
    outdir                  //  string: The output directory where the results will be saved
    input                   //  string: Path to input assemblysheet
    rna_evidence            //  string: Path to rna samplesheet
    liftoff_annotations     //  string: Path to liftoff annotations sheet
    orthofinder_annotations //  string: Path to orthofinder annotations sheet
    help                    // boolean: Display help message and exit
    help_full               // boolean: Show the full help message
    show_hidden             // boolean: Show hidden parameters in the help message

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        true, // validate params
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create input channels
    //
    ch_input                    = Channel.fromList (samplesheetToList(input, "assets/schema_input.json"))

    ch_target_assembly          = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def fasta       = it[1]

                                    def fasta_file  = file(fasta, checkIfExists: true)

                                    if ( workflow.stubRun ) {
                                        return [ [ id: tag ], fasta_file ]
                                    }

                                    def is_zipped   = fasta.endsWith('.gz')
                                    def sz_thresh   = is_zipped ? 300_000 : 1_000_000
                                    def fasta_size  = fasta_file.size()

                                    if ( fasta_size < sz_thresh ) { // < 1 MB
                                        error "The assembly represented by tag '$tag' is only $fasta_size bytes. The minimum allowed size is 1 MB!"
                                    }

                                    [ [ id: tag ], fasta_file ]
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

    ch_is_masked                = ch_input
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

    ch_benchmark_gff            = ch_input
                                | map { it ->
                                    def tag         = it[0]
                                    def gff         = it[6]

                                    if ( gff ) {
                                        [
                                            [ id: tag ],
                                            file(gff, checkIfExists: true)
                                        ]
                                    }
                                }

    ch_rna_branch               = ! params.rna_evidence
                                ? Channel.empty()
                                : Channel.fromList (samplesheetToList(rna_evidence, "assets/schema_rna.json"))
                                | map { meta, f1, f2 ->
                                    "$f1".find(/^SRR[0-9]*$/)
                                    ? [ meta + [ single_end: false ], [ f1, f2 ] ]
                                    : f2
                                    ? [ meta + [ single_end: false ], [ f1, f2 ] ]
                                    : [ meta + [ single_end: true ], [ f1 ] ]
                                }
                                | map { meta, data ->
                                    [ meta + [ target_assemblies: meta.target_assemblies.split(';').sort() ], data ]
                                }
                                | branch { meta, data ->
                                    sra: "${data.first()}".find(/^SRR[0-9]*$/)
                                    fq:  "${data.first()}".find(/^\S+\.(f(ast)?q\.gz)$/)
                                    bam: "${data.first()}".find(/^\S+\.bam$/)
                                }

    ch_rna_sra                  = ! params.rna_evidence
                                ? Channel.empty()
                                : ch_rna_branch.sra
                                | map { meta, sra -> [ meta.id, meta, sra ] }
                                | groupTuple
                                | combine(ch_tar_assm_str)
                                | map { id, metas, sra, tar_assm_str ->
                                    validateSRAMetadata(metas, sra, tar_assm_str)
                                }

    ch_rna_fq                   = ! params.rna_evidence
                                ? Channel.empty()
                                : ch_rna_branch.fq
                                | map { meta, files -> [ meta.id, meta, files.collect { file(it) } ] }
                                | groupTuple
                                | combine(ch_tar_assm_str)
                                | map { id, metas, files, tar_assm_str ->
                                    validateFastqMetadata(metas, files, tar_assm_str)
                                }

    ch_rna_bam                  = ! params.rna_evidence
                                ? Channel.empty()
                                : ch_rna_branch.bam
                                | map { meta, files -> [ meta.id, meta, files.collect { file(it) } ] }
                                | groupTuple
                                | combine(ch_tar_assm_str)
                                | flatMap { id, metas, files, tar_assm_str ->
                                    validateBamMetadata(metas, files, tar_assm_str)
                                }

    // Check if each sample for a given assembly has either bam or fastq files
    ch_rna_bam
    | flatMap { meta, bams ->
        meta.target_assemblies.collect { [ [ meta.id, it ], 'bam' ] }
    }
    | join(
        ch_rna_fq
        | flatMap { meta, fqs ->
            meta.target_assemblies.collect { [ [ meta.id, it ], 'fq' ] }
        }
    )
    | map { combination, bam, fq ->
        error "Sample ${combination[0]} for assembly ${combination[1]} can not have both fastq and bam files"
    }

    ch_rna_bam_by_assembly      = ch_rna_bam
                                | map { meta, bams -> [ [ id: meta.target_assemblies.first() ], bams ] }
                                | groupTuple
                                | map { meta, bams -> [ meta, bams.flatten() ] }

    ch_ribo_db                  = params.remove_ribo_rna
                                ? file(params.ribo_database_manifest, checkIfExists: true)
                                : null

    ch_sortmerna_fastas         = ch_ribo_db
                                ? Channel.from(ch_ribo_db ? ch_ribo_db.readLines() : null)
                                | map { row -> file(row, checkIfExists: true) }
                                | collect
                                : Channel.empty()

    ch_ext_prot_fastas          = ( params.protein_evidence.endsWith('txt')
                                    ? Channel.fromPath(params.protein_evidence)
                                    | splitText
                                    : Channel.fromPath(params.protein_evidence)
                                )
                                | map { file_path ->

                                    def file_handle = ( file_path instanceof String )
                                        ? file(file_path.strip(), checkIfExists: true)
                                        : file_path

                                    [ [ id: idFromFileName( file_handle.baseName ) ], file_handle ]
                                }


    ch_liftoff_mm               = ! params.liftoff_annotations
                                ? Channel.empty()
                                : Channel.fromList (samplesheetToList(liftoff_annotations, "assets/schema_liftoff.json"))
                                | multiMap { fasta, gff ->
                                    def fastaFile = file(fasta, checkIfExists:true)

                                    fasta: [ [ id: idFromFileName( fastaFile.baseName ) ], fastaFile ]
                                    gff: [ [ id: idFromFileName( fastaFile.baseName ) ], file(gff, checkIfExists:true) ]
                                }

    ch_liftoff_fasta            = params.liftoff_annotations
                                ? ch_liftoff_mm.fasta
                                : Channel.empty()

    ch_liftoff_gff              = params.liftoff_annotations
                                ? ch_liftoff_mm.gff
                                : Channel.empty()

    ch_tsebra_config            = Channel.of ( file("${projectDir}/assets/tsebra-template.cfg", checkIfExists: true) )
                                | map { cfg ->
                                    def param_intron_support = params.enforce_full_intron_support ? '1.0' : '0.0'

                                    def param_e1 = params.allow_isoforms ? '0.1'    : '0.0'
                                    def param_e2 = params.allow_isoforms ? '0.5'    : '0.0'
                                    def param_e3 = params.allow_isoforms ? '0.05'   : '0.0'
                                    def param_e4 = params.allow_isoforms ? '0.2'    : '0.0'

                                    [
                                        'tsebra-config.cfg',
                                        cfg
                                        .text
                                        .replace('PARAM_INTRON_SUPPORT', param_intron_support)
                                        .replace('PARAM_E1', param_e1)
                                        .replace('PARAM_E2', param_e2)
                                        .replace('PARAM_E3', param_e3)
                                        .replace('PARAM_E4', param_e4)
                                    ]
                                }
                                | collectFile


    ch_orthofinder_pep          = ! params.orthofinder_annotations
                                ? Channel.empty()
                                : Channel.fromList (samplesheetToList(orthofinder_annotations, "assets/schema_orthofinder.json"))
                                | map { tag, fasta ->
                                    [ [ id: tag ], file(fasta, checkIfExists:true)  ]
                                }

    emit:
    target_assembly             = ch_target_assembly
    is_masked                   = ch_is_masked
    te_library                  = ch_te_library
    braker_annotation           = ch_braker_annotation
    braker_ex_asm_str           = ch_braker_ex_asm_str
    benchmark_gff               = ch_benchmark_gff
    rna_sra                     = ch_rna_sra
    rna_fq                      = ch_rna_fq
    rna_bam_by_assembly         = ch_rna_bam_by_assembly
    sortmerna_fastas            = ch_sortmerna_fastas
    ext_prot_fastas             = ch_ext_prot_fastas
    liftoff_fasta               = ch_liftoff_fasta
    liftoff_gff                 = ch_liftoff_gff
    tsebra_config               = ch_tsebra_config
    orthofinder_pep             = ch_orthofinder_pep
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Additional validation
//
def idFromFileName(fileName) {

    def trial = ( fileName
        ).replaceFirst(
            /\.f(ast)?q$/, ''
        ).replaceFirst(
            /\.f(asta|sa|a|as|aa|na)?$/, ''
        ).replaceFirst(
            /\.gff(3)?$/, ''
        ).replaceFirst(
            /\.gz$/, ''
        )

    if ( trial == fileName ) { return fileName }

    return idFromFileName ( trial )
}

def validateSRAMetadata(metas, sra, permAssString) {
    def permAssList = permAssString.split(",")

    // Check if a single sample is accociated with each SRA ID
    if ( metas.size() > 1 ) {
        error "Samples with SRA IDs should be unique. ${metas.first().id} is repeated"
    }

    def meta = metas.first()

    // Check if sra list has only one element
    if ( sra.first()[1] != [] ) {
        error "Sample '${meta.id}' lists a SRA ID under file1. It's file2 column should be empty"
    }

    // Check if each listed assembly is permissible
    if ( meta.target_assemblies.any { !permAssList.contains( it ) } ) {
        error "Sample ${meta.id} targets ${meta.target_assemblies} which are not in $permAssList"
    }

    [ metas.first(), sra.first().first() ]
}

def validateFastqMetadata(metas, fqs, permAssString) {
    def permAssList = permAssString.split(",")

    // Check if each listed assembly is permissible
    metas.each { meta ->
        if ( meta.target_assemblies.any { !permAssList.contains( it ) } ) {
            error "Sample ${meta.id} targets ${meta.target_assemblies} which are not in $permAssList"
        }
    }

    // Check if multiple runs of a sample have the same target assemblies
    if ( metas.collect { meta -> meta.target_assemblies }.unique().size() > 1 ) {
        error "Multiple runs of sample ${metas.first().id} must target same assemblies"
    }

    // Check if multiple runs of a sample have the same endedness
    if ( metas.collect { meta -> meta.single_end }.unique().size() > 1 ) {
        error "Multiple runs of sample ${metas.first().id} must have same endedness"
    }

    [ metas.first(), fqs ]
}


def validateBamMetadata(metas, bams, permAssString) {
    def permAssList = permAssString.split(",")

    // Check if each listed assembly is permissible
    metas.each { meta ->
        if ( meta.target_assemblies.any { !permAssList.contains( it ) } ) {
            error "Sample ${meta.id} targets ${meta.target_assemblies} which are not in $permAssList"
        }
    }

    // Check that when the first file is bam then the second file is absent
    bams.findAll { files ->
        files.first().extension == 'bam' && files.size() != 1
    }
    .each { error "Sample ${metas.first().id} contains both bam and fastq pairs. When a bam file is provided as file_1, a fastq for file_2 is not permitted" }

    // Check that a bam file only targets a single assembly
    bams.eachWithIndex { files, index ->
        if ( files.first().extension == 'bam' && metas[index].target_assemblies.size() > 1 ) {
            error "BAM file for sample ${metas.first().id} can only target one assembly: ${metas[index].target_assemblies}"
        }
    }

    metas.every { it.target_assemblies == metas.first().target_assemblies }
    ? [ [ metas.first(), bams.flatten() ] ]
    : metas.withIndex().collect { meta, index -> [ meta, bams[index].flatten() ] }
}


// Generate methods description for MultiQC
//
def toolCitationText(versions_yml) {

    def v_text      = versions_yml.text.toLowerCase()

    def start_text  = 'Tools used in the workflow included: '
    def end_text    = ' and MultiQC (Ewels et al. 2016).'

    def citation_text = [
            false                                   ? ''    : 'AGAT (Dainat et al. 2024)',
            false                                   ? ''    : 'AUGUSTUS (Sommerfeld et al. 2009)',
            false                                   ? ''    : 'BRAKER3 (Gabriel et al. 2023)',
            ( ! v_text.contains('busco:') )         ? ''    : 'BUSCO (Manni et al. 2021)',
            ( ! v_text.contains('edta:') )          ? ''    : 'EDTA (Ou et al. 2019)',
            ( ! v_text.contains('eggnog-mapper:') ) ? ''    : 'EggNOG-mapper (Carlos et al. 2021)',
            ( ! v_text.contains('fastqc:') )        ? ''    : 'FASTQC (Andrews. 2010)',
            ( ! v_text.contains('fastp:') )         ? ''    : 'FASTP (Chen et al. 2018)',
            false                                   ? ''    : 'GeneMark-ETP (Brůna et al. 2024)',
            false                                   ? ''    : 'GenomeTools (Gremme et al. 2013)',
            ( ! v_text.contains('gffcompare:') )    ? ''    : 'GffCompare (Pertea & Pertea. 2020)',
            false                                   ? ''    : 'GffRead (Pertea & Pertea. 2020)',
            ( ! v_text.contains('liftoff:') )       ? ''    : 'Liftoff (Shumate & Salzberg. 2021)',
            ( ! v_text.contains('orthofinder:') )   ? ''    : 'OrthoFinder (Emms & Kelly. 2019)',
            false                                   ? ''    : 'ProtHint (Brůna et al. 2020)',
            false                                   ? ''    : 'py_fasta_validator (Edwards. 2019)',
            ( ! v_text.contains('repeatmasker:') )  ? ''    : 'RepeatMasker (Smit & Hubley. 2023)',
            ( ! v_text.contains('repeatmodeler:') ) ? ''    : 'RepeatModeler (Hubley. 2023)',
            false                                   ? ''    : 'Samtools (Danecek et al. 2021)',
            false                                   ? ''    : 'SeqKit (Shen et al. 2016)',
            ( ! v_text.contains('sortmerna:') )     ? ''    : 'SortMeRNA (Kopylova et al. 2012)',
            ( ! v_text.contains('star:') )          ? ''    : 'STAR (Dobin et al. 2013)',
            false                                   ? ''    : 'TSEBRA (Gabriel et al. 2021)',
        ].findAll { it != '' }.join(', ').trim()

    return start_text + citation_text + end_text
}

def toolBibliographyText(versions_yml) {

    def v_text = versions_yml.text.toLowerCase()

    def reference_text = [
            false                                   ? ''    : 'Jacques Dainat, Darío Hereñú, Dr. K. D. Murray, Ed Davis, Ivan Ugrin, Kathryn Crouch, LucileSol, Nuno Agostinho, pascal-git, Zachary Zollman, & tayyrov. (2024). NBISweden/AGAT: AGAT-v1.4.1 (v1.4.1). Zenodo. <a href="https://doi.org/10.5281/zenodo.13799920">10.5281/zenodo.13799920</a>',
            false                                   ? ''    : 'Sommerfeld, D., Lingner, T., Stanke, M., Morgenstern, B., & Richter, H. (2009). AUGUSTUS at MediGRID: Adaption of a bioinformatics application to grid computing for efficient genome analysis. Future Gener. Comput. Syst., 25, 337-345. doi: <a href="https://doi.org/10.1016/j.future.2008.05.010">10.1016/j.future.2008.05.010</a>',
            false                                   ? ''    : 'Gabriel, L., Bruna, T., Hoff, K. J., Ebel, M., Lomsadze, A., Borodovsky, M., Stanke, M. (2023). BRAKER3: Fully Automated Genome Annotation Using RNA-Seq and Protein Evidence with GeneMark-ETP, AUGUSTUS and TSEBRA. bioRxiV, doi: <a href="https://doi.org/10.1101/2023.06.10.544449">10.1101/2023.06.10.544449</a>',
            ( ! v_text.contains('busco:') )         ? ''    : 'Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. 2021. BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes, Molecular Biology and Evolution, Volume 38, Issue 10, October 2021, Pages 4647–4654, doi: <a href="https://doi.org/10.1093/molbev/msab199">10.1093/molbev/msab199</a>',
            ( ! v_text.contains('edta:') )          ? ''    : 'Ou S., Su W., Liao Y., Chougule K., Agda J. R. A., Hellinga A. J., Lugo C. S. B., Elliott T. A., Ware D., Peterson T., Jiang N., Hirsch C. N. and Hufford M. B. (2019). Benchmarking Transposable Element Annotation Methods for Creation of a Streamlined, Comprehensive Pipeline. Genome Biol. 20(1): 275. doi: <a href="https://doi.org/10.1186/s13059-019-1905-y">10.1186/s13059-019-1905-y</a>',
            ( ! v_text.contains('eggnog-mapper:') ) ? ''    : 'eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. Carlos P. Cantalapiedra, Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021. Molecular Biology and Evolution, msab293, doi: <a href="https://doi.org/10.1093/molbev/msab293">10.1093/molbev/msab293</a>',
            ( ! v_text.contains('fastqc:') )        ? ''    : 'Andrews, S. (2010). Babraham Bioinformatics - FastQC a quality control tool for high throughput sequence data. Babraham.ac.uk. url: <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc">https://www.bioinformatics.babraham.ac.uk/projects/fastqc</a>',
            ( ! v_text.contains('fastp:') )         ? ''    : 'Chen S, Zhou Y, Chen Y, Gu J. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 01 September 2018, Pages i884–i890, doi: <a href="https://doi.org/10.1093/bioinformatics/bty560">10.1093/bioinformatics/bty560</a>',
            false                                   ? ''    : 'Brůna T, Lomsadze A, Borodovsky M. GeneMark-ETP significantly improves the accuracy of automatic annotation of large eukaryotic genomes. Genome Res. 2024 Jun 25;34(5):757-768. doi: <a href="https://doi.org/10.1101/gr.278373.123">10.1101/gr.278373.123</a>. PMID: 38866548; PMCID: PMC11216313.',
            false                                   ? ''    : 'Gremme G, Steinbiss S, Kurtz S. 2013. "GenomeTools: A Comprehensive Software Library for Efficient Processing of Structured Genome Annotations," in IEEE/ACM Transactions on Computational Biology and Bioinformatics, vol. 10, no. 3, pp. 645-656, May 2013, doi: <a href="https://doi.org/10.1109/TCBB.2013.68">10.1109/TCBB.2013.68</a>',
            ( ! v_text.contains('gffcompare:') )    ? ''    : 'Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare. F1000Res. 2020 Apr 28;9:ISCB Comm J-304. doi: <a href="http://doi.org/10.12688/f1000research.23297.2">10.12688/f1000research.23297.2</a>. PMID: 32489650; PMCID: PMC7222033.',
            ( ! v_text.contains('liftoff:') )       ? ''    : 'Shumate A, Salzberg SL. Liftoff: accurate mapping of gene annotations. Bioinformatics. 2021 Jul 19;37(12):1639-1643. doi: <a href="http://doi.org/10.1093/bioinformatics/btaa1016">10.1093/bioinformatics/btaa1016</a>. PMID: 33320174; PMCID: PMC8289374.',
            ( ! v_text.contains('orthofinder:') )   ? ''    : 'Emms, D.M., Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol 20, 238 (2019). doi: <a href="https://doi.org/10.1186/s13059-019-1832-y">10.1186/s13059-019-1832-y</a>',
            false                                   ? ''    : 'Tomáš Brůna, Alexandre Lomsadze, Mark Borodovsky, GeneMark-EP+: eukaryotic gene prediction with self-training in the space of genes and proteins, NAR Genomics and Bioinformatics, Volume 2, Issue 2, June 2020, lqaa026, doi: <a href="https://doi.org/10.1093/nargab/lqaa026">10.1093/nargab/lqaa026</a>',
            false                                   ? ''    : 'Edwards, R.A. 2019. fasta_validate: a fast and efficient fasta validator written in pure C. doi: <a href="https://doi.org/10.5281/zenodo.2532044">10.5281/zenodo.2532044</a>',
            ( ! v_text.contains('repeatmasker:') )  ? ''    : 'Smit, A., & Hubley, R. (2023). RepeatMasker: a program that screens DNA sequences for interspersed repeats and low complexity DNA sequences. Repeatmasker.org. url: <a href="https://www.repeatmasker.org">https://www.repeatmasker.org</a>',
            ( ! v_text.contains('repeatmodeler:') ) ? ''    : 'Hubley, R. (2023). RepeatModeler: a de novo transposable element (TE) family identification and modeling package. Repeatmasker.org. url: <a href="https://www.repeatmasker.org">https://www.repeatmasker.org</a>',
            false                                   ? ''    : 'Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. 2021. Twelve years of SAMtools and BCFtools, GigaScience, Volume 10, Issue 2, February 2021, giab008, doi: <a href="https://doi.org/10.1093/gigascience/giab008">10.1093/gigascience/giab008</a>',
            false                                   ? ''    : 'Shen W, Le S, Li Y, Hu F. 2016. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS ONE 11(10): e0163962. doi: <a href="https://doi.org/10.1371/journal.pone.0163962">10.1371/journal.pone.0163962</a>',
            ( ! v_text.contains('sortmerna:') )     ? ''    : 'Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: <a href="https://doi.org/10.1093/bioinformatics/bts611">10.1093/bioinformatics/bts611</a>',
            ( ! v_text.contains('star:') )          ? ''    : 'Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: <a href="https://doi.org/10.1093/bioinformatics/bts635">10.1093/bioinformatics/bts635</a>. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.',
            false                                   ? ''    : 'Gabriel, L., Hoff, K.J., Brůna, T. et al. TSEBRA: transcript selector for BRAKER. BMC Bioinformatics 22, 566 (2021). doi: <a href="https://doi.org/10.1186/s12859-021-04482-0">10.1186/s12859-021-04482-0</a>',
            false                                   ? ''    : 'Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a>',
        ].findAll { it != '' }.collect { it -> "<li>$it</li>" }.join('').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml, versions_yml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText(versions_yml).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText(versions_yml)


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
