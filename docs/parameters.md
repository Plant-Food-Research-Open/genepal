# plantandfoodresearch/pangene pipeline parameters

A NextFlow pipeline for pan-genome annotation

## Input/output options

| Parameter                 | Description                                                                                                                                                                               | Type      | Default   | Required | Hidden |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | --------- | -------- | ------ |
| `input`                   | Target assemblies listed in a CSV sheet <details><summary>Help</summary><small>FASTA and other associated files for target assemblies provided as a CSV sheet</small></details>           | `string`  |           | True     |        |
| `external_protein_fastas` | External protein fastas listed in a text sheet <details><summary>Help</summary><small>A text file listing FASTA files to provide protein evidence for annotation</small></details>        | `string`  |           | True     |        |
| `eggnogmapper_db_dir`     | Eggnogmapper database directory                                                                                                                                                           | `string`  |           | True     |        |
| `eggnogmapper_tax_scope`  | Eggnogmapper taxonomy scopre                                                                                                                                                              | `integer` |           | True     |        |
| `fastq`                   | FASTQ samples listed in a CSV sheet <details><summary>Help</summary><small>FASTQ files for RNASeq samples corresponding to each target assembly provided in a CSV sheet</small></details> | `string`  |           |          |        |
| `liftoff_annotations`     | Reference annotations listed in a CSV sheet <details><summary>Help</summary><small>FASTA and GFF3 files for reference annotations for liftoff listed in a CSV sheet</small></details>     | `string`  |           |          |        |
| `outdir`                  | The output directory where the results will be saved <details><summary>Help</summary><small> Use absolute paths to storage on Cloud infrastructure</small></details>                      | `string`  | ./results | True     |        |

## Repeat annotation options

| Parameter                   | Description                                | Type      | Default       | Required | Hidden |
| --------------------------- | ------------------------------------------ | --------- | ------------- | -------- | ------ |
| `repeat_annotator`          | 'edta' or 'repeatmodeler'                  | `string`  | repeatmodeler |          |        |
| `save_annotated_te_lib`     | Save annotated TE library or not?          | `boolean` |               |          |        |
| `edta_is_sensitive`         | Use '--sensitive 1' flag with EDTA or not? | `boolean` |               |          |        |
| `repeatmasker_save_outputs` | Save the repeat-masked genome or not?      | `boolean` |               |          |        |

## RNASeq pre-processing options

| Parameter                | Description                                                        | Type      | Default                                   | Required | Hidden |
| ------------------------ | ------------------------------------------------------------------ | --------- | ----------------------------------------- | -------- | ------ |
| `skip_fastqc`            | Skip FASTQC or not?                                                | `boolean` |                                           |          |        |
| `skip_fastp`             | Skip trimming by FASTQP or not?                                    | `boolean` |                                           |          |        |
| `min_trimmed_reads`      | Exclude a sample if its reads after trimming are below this number | `integer` | 10000                                     |          |        |
| `extra_fastp_args`       | Extra FASTP arguments                                              | `string`  |                                           |          |        |
| `save_trimmed`           | Save FASTQ files after trimming or not?                            | `boolean` |                                           |          |        |
| `remove_ribo_rna`        | Remove Ribosomal RNA or not?                                       | `boolean` |                                           |          |        |
| `save_non_ribo_reads`    | Save FASTQ files after Ribosomal RNA removal or not?               | `boolean` |                                           |          |        |
| `ribo_database_manifest` | Ribosomal RNA fastas listed in a text sheet                        | `string`  | ${projectDir}/assets/rrna-db-defaults.txt |          |        |

## RNAseq alignment options

| Parameter                | Description                                       | Type      | Default | Required | Hidden |
| ------------------------ | ------------------------------------------------- | --------- | ------- | -------- | ------ |
| `star_max_intron_length` | Maximum intron length for STAR alignment          | `integer` | 16000   |          |        |
| `star_align_extra_args`  | EXTRA arguments for STAR                          | `string`  |         |          |        |
| `star_save_outputs`      | Save BAM files from STAR or not?                  | `boolean` |         |          |        |
| `save_cat_bam`           | SAVE a concatenated BAM file per assembly or not? | `boolean` |         |          |        |

## Annotation options

| Parameter                   | Description                                                                       | Type      | Default | Required | Hidden |
| --------------------------- | --------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `braker_extra_args`         | Extra arguments for BRAKER                                                        | `string`  |         |          |        |
| `braker_allow_isoforms`     | Allow multiple isoforms for gene models                                           | `boolean` | True    |          |        |
| `liftoff_coverage`          | Liftoff coverage parameter                                                        | `number`  | 0.9     |          |        |
| `liftoff_identity`          | Liftoff identity parameter                                                        | `number`  | 0.9     |          |        |
| `eggnogmapper_evalue`       | Only report alignments below or equal the e-value threshold                       | `number`  | 1e-05   |          |        |
| `eggnogmapper_pident`       | Only report alignments above or equal to the given percentage of identity (0-100) | `integer` | 35      |          |        |
| `eggnogmapper_purge_nohits` | Purge transcripts which do not have a hit against eggnog                          | `boolean` |         |          |        |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter    | Description                                                                                                                                                                                                                                                                 | Type      | Default | Required | Hidden |
| ------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `max_cpus`   | Maximum number of CPUs that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the CPU requirement for each process. Should be an integer e.g. `--max_cpus 1`</small></details>                                      | `integer` | 12      |          | True   |
| `max_memory` | Maximum amount of memory that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the memory requirement for each process. Should be a string in the format integer-unit e.g. `--max_memory '8.GB'`</small></details> | `string`  | 200.GB  |          | True   |
| `max_time`   | Maximum amount of time that can be requested for any single job. <details><summary>Help</summary><small>Use to set an upper-limit for the time requirement for each process. Should be a string in the format integer-unit e.g. `--max_time '2.h'`</small></details>        | `string`  | 7.day   |          | True   |
