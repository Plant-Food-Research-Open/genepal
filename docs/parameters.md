

# plant-food-research-open/genepal pipeline parameters

A Nextflow pipeline for consensus, phased and pan-genome annotation.

## Input/output options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Target assemblies listed in a CSV sheet | `string` |  | True |  |
| `protein_evidence` | Protein evidence provided as a fasta file or multiple fasta files listed in a plain txt file | `string` |  | True |  |
| `eggnogmapper_db_dir` | Eggnogmapper database directory | `string` |  |  |  |
| `eggnogmapper_tax_scope` | Eggnogmapper taxonomy scopre. Eukaryota: 2759, Viridiplantae: 33090, Archaea: 2157, Bacteria: 2, root: 1 | `integer` | 1 |  |  |
| `min_contig_length` | Minimum length for contigs to be included in annotation | `integer` | 5000 |  |  |
| `rna_evidence` | FASTQ/BAM samples listed in a CSV sheet | `string` |  |  |  |
| `liftoff_annotations` | Reference annotations listed in a CSV sheet | `string` |  |  |  |
| `orthofinder_annotations` | Additional annotations for orthology listed in a CSV sheet | `string` |  |  |  |
| `outdir` | The output directory where the results will be saved | `string` |  | True |  |
| `email` | Email address for completion summary. | `string` |  |  | True |

## Repeat annotation options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `repeat_annotator` | 'edta' or 'repeatmodeler' | `string` | repeatmodeler |  |  |
| `save_annotated_te_lib` | Save annotated TE library or not? | `boolean` |  |  |  |
| `edta_is_sensitive` | Use '--sensitive 1' flag with EDTA or not? | `boolean` |  |  |  |
| `repeatmasker_save_outputs` | Save the repeat-masked genome or not? | `boolean` |  |  |  |

## RNASeq pre-processing options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fastqc_skip` | Skip FASTQC or not? | `boolean` | True |  |  |
| `fastp_skip` | Skip trimming by FASTQP or not? | `boolean` |  |  |  |
| `min_trimmed_reads` | Exclude a sample if its reads after trimming are below this number | `integer` | 10000 |  |  |
| `fastp_extra_args` | Extra FASTP arguments | `string` |  |  |  |
| `save_trimmed` | Save FASTQ files after trimming or not? | `boolean` |  |  |  |
| `remove_ribo_rna` | Remove Ribosomal RNA or not? | `boolean` |  |  |  |
| `save_non_ribo_reads` | Save FASTQ files after Ribosomal RNA removal or not? | `boolean` |  |  |  |
| `ribo_database_manifest` | Ribosomal RNA fastas listed in a text sheet | `string` | ${projectDir}/assets/rrna-db-defaults.txt |  |  |

## RNASeq alignment options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `star_max_intron_length` | Maximum intron length for STAR alignment | `integer` | 16000 |  |  |
| `star_align_extra_args` | EXTRA arguments for STAR | `string` |  |  |  |
| `star_save_outputs` | Save BAM files from STAR or not? | `boolean` |  |  |  |
| `save_cat_bam` | SAVE a concatenated BAM file per assembly or not? | `boolean` |  |  |  |

## Annotation options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `braker_extra_args` | Extra arguments for BRAKER | `string` |  |  |  |
| `liftoff_coverage` | Liftoff coverage parameter | `number` | 0.9 |  |  |
| `liftoff_identity` | Liftoff identity parameter | `number` | 0.9 |  |  |
| `eggnogmapper_evalue` | Only report alignments below or equal the e-value threshold | `number` | 1e-05 |  |  |
| `eggnogmapper_pident` | Only report alignments above or equal to the given percentage of identity (0-100) | `integer` | 35 |  |  |

| `--min_contig_length` | Minimum length (in base pairs) of contigs to include in the analysis    | `number`| 5000    |          |        |


## Post-annotation filtering options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `allow_isoforms` | Allow multiple isoforms for gene models | `boolean` | True |  |  |
| `enforce_full_intron_support` | Require every model to have external evidence for all its introns | `boolean` | True |  |  |
| `filter_liftoff_by_hints` | Use BRAKER hints to filter Liftoff models | `boolean` | True |  |  |
| `eggnogmapper_purge_nohits` | Purge transcripts which do not have a hit against eggnog | `boolean` |  |  |  |

## Annotation output options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `braker_save_outputs` | Save BRAKER files | `boolean` |  |  |  |
| `add_attrs_to_proteins_fasta` | Add gff attributes to proteins fasta | `boolean` |  |  |  |

## Evaluation options



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `busco_skip` | Skip evaluation by BUSCO | `boolean` |  |  |  |
| `busco_lineage_datasets` | BUSCO lineages as a space-separated list: 'fungi_odb10 microsporidia_odb10' | `string` | eukaryota_odb10 |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. | `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. | `string` | copy |  | True |
| `email_on_fail` | Email address for completion summary, only when pipeline fails. | `string` |  |  | True |
| `plaintext_email` | Send plain-text email instead of HTML. | `boolean` |  |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `hook_url` | Incoming hook URL for messaging service | `string` |  |  | True |



