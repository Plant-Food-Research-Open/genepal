# plant-food-research-open/genepal: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

<!-- no toc -->

- [Repeat annotation](#repeat-annotation)
- [Repeat masking](#repeat-masking)
- [RNASeq trimming, filtering and QC](#rnaseq-trimming-filtering-and-qc)
- [RNASeq alignment](#rnaseq-alignment)
- [Annotation with BRAKER](#annotation-with-braker)
- [Annotation with Liftoff](#annotation-with-liftoff)
- [Annotation filtering and merging](#annotation-filtering-and-merging)
- [Functional annotation](#functional-annotation)
- [Orthology inference](#orthology-inference)
- [Final annotation files](#final-annotation-files)
- [Annotation QC](#annotation-qc)
- [Pipeline information and MultiQC](#pipeline-information-and-multiqc)

### Repeat annotation

<details markdown="1">
<summary>Output files</summary>

- `repeatmodeler/`
  - `*.fa`: Repeat library
- `edta/`
  - `*.EDTA.TElib.fa`: Repeat library

</details>

A repeat library is created with either [REPEATMODELER](https://github.com/Dfam-consortium/RepeatModeler) or [EDTA](https://github.com/oushujun/EDTA). The choice of the tool is specified by the `repeat_annotator` parameter (default: `repeatmodeler`). Repeat annotation outputs are saved to the output directory only if `save_annotated_te_lib` parameter is set to `true` (default: `false`).

### Repeat masking

<details markdown="1">
<summary>Output files</summary>

- `repeatmasker/`
  - `*.masked`: Masked assembly

</details>

Soft masking of the repeats is performed with [REPEATMASKER](https://github.com/rmhubley/RepeatMasker) using the repeat library prepared in the previous step. Masking outputs are saved to the output directory only if `repeatmasker_save_outputs` parameter is set to `true` (default: `false`).

### RNASeq trimming, filtering and QC

<details markdown="1">
<summary>Output files</summary>

- `fastqc_raw/`
  - `*.html`: HTML QC report for a sample before trimming
  - `*.zip`: Zipped QC files for a sample before trimming
- `fastqc_trim/`
  - `*.html`: HTML QC report for a sample after trimming
  - `*.zip`: Zipped QC files for a sample after trimming
- `fastp/`
  - `html/`
    - `*.fastp.html`: HTML trimming report for a sample
  - `json/`
    - `*.fastp.json`: Trimming statistics for a sample
  - `log/`
    - `*.fastp.log`: Trimming log for a sample
  - `*_{1,2}.fail.fastq.gz`: Reads which failed trimming
  - `*.paired.fail.fastq.gz`: Pairs of reads which failed trimming
  - `*.merged.fastq.gz`: Reads which passed trimming. For paired reads, reads 1 and 2 are merged into a single file
- `sortmerna/`
  - `*.sortmerna.log`: Filtering log for a sample
  - `*_{1,2}.non_rRNA.fastq.gz`: Filtered reads

</details>

RNASeq reads are trimmed with [FASTP](https://github.com/OpenGene/fastp) and are QC'ed with [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc). Ribosomal reads are filtered out using [SORTMERNA](https://github.com/sortmerna/sortmerna). Trimmed reads are only stored to the output directory if the `save_trimmed` parameter is set to `true` (default: `false`). Reads filtered by [SORTMERNA](https://github.com/sortmerna/sortmerna) are stored to the output directory if the `save_non_ribo_reads` parameter is set to `true` (default: `false`).

### RNASeq alignment

<details markdown="1">
<summary>Output files</summary>

- `star/`
  - `alignment/`
    - `X.on.Y.Aligned.sortedByCoord.out.bam`: Sorted BAM file of read alignments for sample `X` against reference `Y`
    - `X.on.Y.Log.final.out`: STAR final log file for sample `X` against reference `Y`
  - `cat_bam/`
    - `Y.bam`: A single BAM file for reference `Y` created by concatenating alignments from sample-wise `*.on.Y.Aligned.sortedByCoord.out.bam` files

</details>

RNASeq alignment is performed with [STAR](https://github.com/alexdobin/STAR). Alignment files are only stored to the output directory if the `star_save_outputs` parameter is set to `true` (default: `false`). Concatenated bam files are stored to the output directory if the `save_cat_bam` parameter is set to `true` (default: `false`).

### Annotation with BRAKER

<details markdown="1">
<summary>Output files</summary>

- `etc/braker/`
  - `Y/`
    - `braker.gff3`: Gene models predicted by BRAKER in GFF3 format
    - `braker.gtf`: Gene models predicted by BRAKER in GTF format
    - `braker.codingseq`: Coding sequences for the predicted genes
    - `braker.aa`: Protein sequences for the predicted genes
    - `braker.log`: BRAKER log file
    - `hintsfile.gff`: Evidential hints used by BRAKER in GFF format
    - `what-to-cite.txt`: A list of references which must be cited when reporting outputs created by BRAKER

</details>

[BRAKER](https://github.com/Gaius-Augustus/BRAKER) is used to annotate each genome assembly using the provide protein and RNASeq evidence. Outputs from BRAKER are stored to the output directory if the `braker_save_outputs` parameter is set to `true` (default: `false`).

> [!CAUTION]
>
> BRAKER outputs are not the final outputs from the pipeline and that's why they are not stored by default. These are only intermediary files.
>
> The pipeline further processes the BRAKER predictions and stores the final validated outputs in the `annotations` directory. The `braker_save_outputs` option is only provided to allow a manual resume of the pipeline for advanced use cases. See [Advanced inputs for manual resume](./usage.md#advanced-inputs-for-manual-resume) in the [usage doc](./usage.md).

### Annotation with Liftoff

Gene models are lifted from reference assembly(ies) to the target assembly using [LIFTOFF](https://github.com/agshumate/Liftoff). Currently, the outputs from Liftoff are considered intermediary and an option to store them in the output directory is not available.

### Annotation filtering and merging

Annotations obtained from [BRAKER](https://github.com/Gaius-Augustus/BRAKER) and [LIFTOFF](https://github.com/agshumate/Liftoff) are filtered with [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA) and merged with [AGAT](https://github.com/NBISweden/AGAT). Currently, the outputs from these processes are considered intermediary and an option to store them in the output directory is not available.

### Functional annotation

<details markdown="1">
<summary>Output files</summary>

- `annotations/`
  - `Y/`
    - `Y.emapper.annotations`: TSV with the annotation results
    - `Y.emapper.hits`: TSV with the search results
    - `Y.emapper.seed_orthologs`: TSV with the results from parsing the hits, linking queries with seed orthologs

</details>

Functional annotation of the gene models from BRAKER and Liftoff is performed with [EGGNOG-MAPPER](https://github.com/eggnogdb/eggnog-mapper).

### Orthology inference

<details markdown="1">
<summary>Output files</summary>

- `orthofinder/`
  - `Comparative_Genomics_Statistics/`
  - `Gene_Duplication_Events/`
  - `Orthogroups/`
  - `Phylogenetic_Hierarchical_Orthogroups/`
  - `Species_Tree/`

</details>

If more than one genome is included in the pipeline, [ORTHOFINDER](https://github.com/davidemms/OrthoFinder) is used to perform an orthology inference.

### Final annotation files

<details markdown="1">
<summary>Output files</summary>

- `annotations/`
  - `Y/`
    - `Y.gt.gff3`: Final annotation file for genome `Y` which contains gene models and their functional annotations
    - `Y.pep.fasta`: Protein sequences for the gene models
    - `Y.cdna.fasta`: cDNA sequences for the gene models
    - `Y.cds.fasta`: Coding sequences for the gene models

</details>

The final annotation files are saved in GFF3 format validated with [GENOMETOOLS](https://github.com/genometools/genometools) and FASTA format obtained with [GFFREAD](https://github.com/gpertea/gffread).

### Annotation QC

<details markdown="1">
<summary>Output files</summary>

- `busco/`
  - `gff/`
    - `short_summary.specific.Y.eudicots_odb10.txt`: BUSCO summary for annotations from genome `Y` against the `eudicots_odb10` database
    - `busco_figure`: BUSCO summary figure including statistics for annotations from all the genomes
  - `fasta/`
    - `short_summary.specific.Y.eudicots_odb10.txt`: BUSCO summary for genome `Y` against the `eudicots_odb10` database
    - `busco_figure`: BUSCO summary figure including statistics for all the genomes
- `etc/`
  - `splicing_marked/`
    - `Y.gff3`: Final annotation file for genome `Y` which contains gene models and their functional annotations. Additionally, the intron features are marked as canonical or non-canonical and the splice motif is also added an attribute.

</details>

The completeness of the annotations is checked with [BUSCO](https://gitlab.com/ezlab/busco). To provide a comparative baseline, the completeness of the genomes is also checked. Moreover, the canonical/non-canonical splicing of the introns is also assessed by the pipeline.

### Reports

<details markdown="1">
<summary>Output files</summary>

- `multiqc_report.html`: A MultiQC report which includes QC statistics, software versions and references
- `genepal_report.html`: A specialized pangene analysis report
- `genepal_data/`: Files containing parsed data from the reporting module

</details>

[MultiQC](https://seqera.io/multiqc/) and [R Markdown](https://rmarkdown.rstudio.com) are used to curate the results of the pipeline in two HTML reports. The MultiQC report is meant to serve as an exhaustive report of the pipeline outputs. Whereas, the R Markdown report (`genepal_report.html`) provides an overall summary along with a pangene analysis if multiple genomes are provided as input.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `genepal_software_mqc_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage. [MultiQC](https://multiqc.info) compiles a HTML report from the tools used by the pipeline.
