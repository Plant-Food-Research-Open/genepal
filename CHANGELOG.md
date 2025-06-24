# plant-food-research-open/genepal: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [24-Jun-2025]

### `Added`

1. Updated the nf-core template to `3.3.1`

### `Dependencies`

1. Nextflow!>=24.04.2
2. nf-schema@2.2.0

## v0.7.2 - [06-May-2025]

### `Fixed`

1. Now using `agat_sp_complement_annotations.pl` to merge Liftoff and BRAKER models to avoid creation of iso-forms due to overlap of separate genes [#153](https://github.com/Plant-Food-Research-Open/genepal/issues/153)

### `Dependencies`

1. Nextflow!>=24.04.2
2. nf-schema@2.2.0

## v0.7.1 - [02-May-2025]

### `Fixed`

1. Gene models from BRAKER with invalid ORF(s) are now removed [#151](https://github.com/Plant-Food-Research-Open/genepal/issues/151)
2. Demoted nf-schema to 2.2.0 to avoid errors with latest Nextflow versions
3. Fixed a nextflow syntax issue in `conf/modules.config`

### `Dependencies`

1. Nextflow!>=24.04.2
2. nf-schema@2.2.0

## v0.7.0 - [16-April-2025]

### `Added`

1. Added parameter `append_genome_prefix_to_feature_ids` which allows the user to add genome prefixes defined in the assemblysheet to the final Gff/Fasta files [#135](https://github.com/Plant-Food-Research-Open/genepal/issues/135)
2. Updated nf-core template to 3.2.0

### `Fixed`

1. Fixed an issue where `filter_genes_by_aa_length` was not correctly applied when the CDS was shorter than the transcript by replacing `GFFREAD` with `AGAT_SPFILTERBYORFSIZE` [#139](https://github.com/Plant-Food-Research-Open/genepal/issues/139)

### `Dependencies`

1. Nextflow!>=24.04.2
2. nf-schema@2.3.0

### `Tool Updates`

| Tool      | Old Version | New Version |
| --------- | ----------- | ----------- |
| agat      | 1.4.0       | 1.4.2       |
| braker3   | v3.0.7.5    | v3.0.7.6    |
| busco     | 5.7.1       | 5.8.3       |
| coreutils | 8.30        | 9.5         |
| fastp     | 0.23.4      | 0.24.0      |
| multiqc   | 1.25.1      | 1.28        |
| seqkit    | 2.8.1       | 2.9.0       |
| htslib    | 1.18        | 1.21        |
| samtools  | 1.18        | 1.21        |
| star      | 2.7.10a     | 2.7.11b     |

## v0.6.0 - [20-Dec-2024]

### 'Added'

1. Added cDNA and CDS outputs to <OUTPUT_DIR>/annotations/<SAMPLE> directory [#118](https://github.com/Plant-Food-Research-Open/genepal/issues/118)
2. Added parameter `add_attrs_to_proteins_cds_fastas`
3. Added parameter `filter_genes_by_aa_length` with default set to `24` which allows removal of genes with ORFs shorter than 24 [#125](https://github.com/Plant-Food-Research-Open/genepal/issues/125)

### `Fixed`

1. Fixed an issue where TSEBRA failed because LIFTOFF lifted non-protein coding genes [#121](https://github.com/Plant-Food-Research-Open/genepal/issues/121)
2. Switched branch name from `master` to `main` in the GHA CIs
3. Fixed an issue in `genepal_report.Rmd` which caused the pangene matrix plot to fail when the number of clusters exceeded 65536 [#124](https://github.com/Plant-Food-Research-Open/genepal/issues/124)
4. Fixed an issue where `GENEPALREPORT` process failed due to OOM kill signal from SLURM [#123](https://github.com/Plant-Food-Research-Open/genepal/issues/123)
5. Fixed an issue where Gff merge after liftoff failed when one of the Gff files did not contain any genes
6. Fixed an issue where `gxf_fasta_agat_spaddintrons_spextractsequences` crashed due to short introns [#89](https://github.com/Plant-Food-Research-Open/genepal/issues/89)

### `Dependencies`

1. Nextflow!>=24.04.2
2. nf-schema@2.1.1

### `Deprecated`

1. Removed parameter `add_attrs_to_proteins_fasta`

## v0.5.0 - [21-Nov-2024]

### `Added`

1. Added MultiQC [#65](https://github.com/plant-food-research-open/genepal/issues/65)
2. Updated nf-core template to 3.0.2 [#66](https://github.com/PlantandFoodResearch/genepal/issues/66)
3. Integrated nf-test into pipeline CI [#68](https://github.com/PlantandFoodResearch/genepal/issues/68)
4. Updated the flowchart [#87](https://github.com/PlantandFoodResearch/genepal/issues/87)
5. Added a large test dataset for the `test_full` profile [#90](https://github.com/PlantandFoodResearch/genepal/issues/90)
6. Now `.gff.gz` and `.gff3.gz` inputs are also allowed for the `benchmark` column in `--input`
7. Now removing liftoff genes with any intron shorted than 10bp [#89](https://github.com/Plant-Food-Research-Open/genepal/issues/89)
8. Now also removing `rRNA` and `tRNA` after liftoff as the downstream logic in the pipeline can not correctly handle these
9. Now skipping FastQC by default [#98](https://github.com/Plant-Food-Research-Open/genepal/issues/98)
10. Added an HTML report [#44](https://github.com/Plant-Food-Research-Open/genepal/issues/44)
11. Added content type as text/html for the MultiQC and genepal reports
12. Added sra-tools for RNASeq data download [#102](https://github.com/Plant-Food-Research-Open/genepal/issues/102)

### `Fixed`

1. Now using `${meta.id}_trim` as prefix for `FASTQC` files
2. Updated citations to include DOIs
3. Fixed a bug where FASTQ versions were not correctly captured
4. Now using the correct out channel from `STAR_ALIGN`. This bug was introduced by a module update during the development of this version [#74](https://github.com/Plant-Food-Research-Open/genepal/issues/74)
5. Fixed OrthoFinder results copy failure on AWS [#108](https://github.com/Plant-Food-Research-Open/genepal/issues/108)

### `Dependencies`

1. Nextflow!>=24.04.2
2. nf-schema@2.1.1

### `Deprecated`

1. Resource parameters have been removed: `max_memory`, `max_cpus`, `max_time`
2. Removed a number of unnecessary parameters: `monochromeLogs`, `config_profile_contact`, `config_profile_url`, `validationFailUnrecognisedParams`, `validationLenientMode`, `validationSchemaIgnoreParams`, `validationShowHiddenParams`, `validate_params`
3. Removed `extra_fastp_args` and replaced it with `fastp_extra_args`
4. Removed and replaced `skip_fastp` and `skip_fastqc` with `fastp_skip` and `fastqc_skip` [#82](https://github.com/Plant-Food-Research-Open/genepal/issues/82)

## v0.4.0 - [04-Oct-2024]

### `Added`

1. Added `orthofinder_annotations` param
2. Added `FASTA_GFF_ORTHOFINDER` sub-workflow
3. Added evaluation by BUSCO [#41](https://github.com/plant-food-research-open/genepal/issues/41)
4. Included common tax ids for eggnog mapper [#27](https://github.com/plant-food-research-open/genepal/issues/27)
5. Implemented hierarchical naming scheme: geneI.tJ, geneI.tJ.exonK, geneI.tJ.cdsK [#19](https://github.com/plant-food-research-open/genepal/issues/19), [#34](https://github.com/plant-food-research-open/genepal/issues/34)
6. Now sorting list of bam and list of fastq before cat to avoid resume cache misses
7. Allowed BAM files for RNA evidence [#3](https://github.com/plant-food-research-open/genepal/issues/3)
8. Added `GXF_FASTA_AGAT_SPADDINTRONS_SPEXTRACTSEQUENCES` sub-workflow for splice type statistics [#11](https://github.com/plant-food-research-open/genepal/issues/11)
9. Changed `orthofinder_annotations` from FASTA/GFF to protein FASTA [#43](https://github.com/plant-food-research-open/genepal/issues/43)
10. Added param `enforce_full_intron_support` to turn on/off strict model purging by TSEBRA [#21](https://github.com/plant-food-research-open/genepal/issues/21)
11. Added param `filter_liftoff_by_hints` to evaluate liftoff models with TSEBRA to make sure they have the same level of evidence as BRAKER [#28](https://github.com/plant-food-research-open/genepal/issues/28)
12. Added a script to automatically check module version updates
13. Reduced `BRAKER3` threads to 8 [#55](https://github.com/plant-food-research-open/genepal/issues/55)
14. Now the final annotations are stored in the `annotations` folder [#53](https://github.com/plant-food-research-open/genepal/issues/53)
15. Now a single `fasta` file can be directly specified for `protein_evidence`
16. `eggnogmapper_db_dir` is not a required parameter anymore
17. `eggnogmapper_tax_scope` is now set to 1 (root div) by default
18. Added a `test` profile based on public data
19. Added parameter `add_attrs_to_proteins_fasta` to enable/disable addition of decoded gff attributes to proteins fasta [#58](https://github.com/plant-food-research-open/genepal/issues/58)
20. Added a check for input assemblies. If an assembly is smaller than 1 MB (or 300KB in zipped format), the pipeline errors out before starting the downstream processes [#47](https://github.com/plant-food-research-open/genepal/issues/47)
21. Now `REPEATMASKER` GFF output is saved via `CUSTOM_RMOUTTOGFF3` [#54](https://github.com/plant-food-research-open/genepal/issues/54)
22. Added `benchmark` column to the input sheet and used `GFFCOMPARE` to perform benchmarking [#63](https://github.com/plant-food-research-open/genepal/issues/63)
23. Added `SEQKIT_RMDUP` to detect duplicate sequence and wrap the fasta to 80 characters
24. Updated parameter section labels for annotation and post-annotation filtering [#64](https://github.com/plant-food-research-open/genepal/issues/64)
25. Updated modules and sub-workflows

### `Fixed`

1. Fixed BRAKER spellings [#36](https://github.com/plant-food-research-open/genepal/issues/36)
2. Fixed liftoff failure when lifting off from a single reference [#40](https://github.com/plant-food-research-open/genepal/issues/40)
3. Added versions from GFF_STORE sub-workflows [#33](https://github.com/plant-food-research-open/genepal/issues/33)

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`

1. Renamed `external_protein_fastas` param to `protein_evidence`
2. Renamed `fastq` param to `rna_evidence`
3. Renamed `braker_allow_isoforms` param to `allow_isoforms`
4. Moved liftoffID from gene level to mRNA/transcript level
5. Moved `version_check.sh` to `.github/version_checks.sh`
6. Removed dependency on <https://github.com/kherronism/nf-modules.git> for `BRAKER3` and `REPEATMASKER` modules which are now installed from <https://github.com/GallVp/nxf-components.git>
7. Removed dependency on <https://github.com/PlantandFoodResearch/nxf-modules.git>
8. Now the final annotations are not stored in the `final` folder
9. Now BRAKER3 outputs are not saved by default [#53](https://github.com/plant-food-research-open/genepal/issues/53) and saved under `etc` folder when enabled
10. Removed `local` profile. Local executor is the default when no executor is specified. Therefore, the `local` profile was not needed.
11. Removed `CUSTOM_DUMPSOFTWAREVERSIONS`
12. `pipeline_info/software_versions.yml` has been replaced with `pipeline_info/genepal_software_mqc_versions.yml`

## v0.3.3 - [18-Jun-2024]

### `Added`

1. Added a stub test to evaluate the case where an assembly is soft masked but has no annotations

### `Fixed`

1. Fixed a bug where `is_masked` was ignored by the pipeline
2. Fixed a bug in param validation which allowed specification of `braker_hints` without `braker_gff3`

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`

## v0.3.2 - [13-May-2024]

### `Added`

### `Fixed`

1. Increased time limit for REPEATMODELER_REPEATMODELER to 5 days
2. Now removing comments from fasta file before feeding it to BRAKER added tests for the perl one liner
3. Fixed CHANGELOG version check failure in `version_check.sh`
4. Increased the SLURM job time limit to 14 days

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`

## v0.3.1 - [10-May-2024]

### `Added`

### `Fixed`

1. Increased time limit for REPEATMODELER_REPEATMODELER to 3 days, REPEATMASKER to 2 days, EDTA_EDTA to 7 days, BRAKER3 to 7 days and EGGNOGMAPPER to 1 day

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`

## v0.3.0 - [30-April-2024]

### `Added`

1. Added changelog and semantic versioning
2. Changed license to MIT
3. Updated `.editorconfig`
4. Moved .literature to test/ branch
5. Renamed `genepal_local` to `local_genepal`
6. Renamed `genepal_pfr` to `pfr_genepal`
7. Added versioning checking
8. Updated github workflow to use pre-commit instead of prettier and editorconfig check
9. Added central singularity cache dir for pfr config
10. Added `SORTMERNA_INDEX` before `SORTMERNA`
11. Fixed sample contamination bug introduced by `file.simpleName`
12. Now using empty files for stub testing in CI
13. Now BRAKER can be skipped by including BRAKER outputs from previous runs in the `target_assemblies` param
14. Added `gffcompare` to merge liftoff annotations
15. Renamed `samplesheet` param to `fastq`
16. Now using assemblysheet in combination with nf-validation for assembly input
17. Added nextflow_schema.json
18. Now using nf-validation to validate fastqsheet provided by params.fastq
19. Moved `manifest.config` and `reporting_defaults.config` content to `nextflow.config`
20. Now using a txt file for `params.external_protein_fastas`
21. Now using nf-validation for `params.liftoff_annotations`
22. Now using nf-validation for all the parameters
23. Added `PURGE_BRAKER_MODELS` sub-workflow
24. Added `GFF_EGGNOGMAPPER` sub-workflow
25. Now using a custom version of `GFFREAD` which supports `meta` and `fasta`
26. Now using TSEBRA to purge models which do not have full intron support from BRAKER hints
27. Added params `eggnogmapper_evalue` and `eggnogmapper_pident`
28. Added `PURGE_NOHIT_BRAKER_MODELS` sub-workflow
29. Now merging BRAKER and liftoff models before running eggnogmapper
30. Added `GFF_MERGE_CLEANUP` sub-workflow
31. Now using `description` field to store notes and textual annotations in the gff files
32. Now using `mRNA` in place of `transcript` in gff files
33. Now `eggnogmapper_purge_nohits` is set to `false` by default
34. Added `GFF_STORE` sub workflow
35. `external_protein_fastas` and `eggnogmapper_db_dir` are not mandatory parameters
36. Added contributors
37. Add a document for the pipeline parameters
38. Updated `pfr_genepal` and `pfr/profile.config`
39. Now using local tests/stub files for GitHub CI
40. Now removing iso-forms left by TSEBRA using `AGAT_SPFILTERFEATUREFROMKILLLIST`
41. Added `pyproject.toml`
42. Now using PFAMs from eggnog if description is '-'

### `Fixed`

1. Removed liftoff models with `valid_ORF=False`
2. Updated license text to include 'Copyright (c) 2024 The New Zealand Institute for Plant and Food Research Limited'

### `Dependencies`

1. NextFlow!>=23.04.4
2. nf-validation=1.1.3

### `Deprecated`
