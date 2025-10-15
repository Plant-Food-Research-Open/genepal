# plant-food-research-open/genepal: Usage<!-- omit in toc -->

> [!NOTE]
>
> This document does not describe every pipeline parameter. For an exhaustive list of parameters, see [parameters.md](./parameters.md).

- [Assemblysheet input](#assemblysheet-input)
  - [Advanced inputs for manual resume](#advanced-inputs-for-manual-resume)
- [Protein evidence](#protein-evidence)
  - [BRAKER workflow](#braker-workflow)
- [RNASeq evidence](#rnaseq-evidence)
  - [BRAKER workflow](#braker-workflow-1)
  - [Preprocessing](#preprocessing)
  - [Alignment](#alignment)
- [Liftoff annotations](#liftoff-annotations)
- [EggNOG-mapper DB](#eggnog-mapper-db)
- [Orthology inference input](#orthology-inference-input)
- [Iso-forms and full intron support](#iso-forms-and-full-intron-support)
- [Running the pipeline](#running-the-pipeline)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Core Nextflow arguments](#core-nextflow-arguments)
  - [`-profile`](#-profile)
  - [`-resume`](#-resume)
  - [`-c`](#-c)
- [Custom configuration](#custom-configuration)
  - [Resource requests](#resource-requests)
  - [Custom Containers](#custom-containers)
  - [Custom Tool Arguments](#custom-tool-arguments)
  - [nf-core/configs](#nf-coreconfigs)
- [Running in the background](#running-in-the-background)
- [Nextflow memory requirements](#nextflow-memory-requirements)

## Assemblysheet input

> ✅ Mandatory `--input`

You will need to create an assemblysheet with information about the genome assemblies you would like to annotate before running the pipeline. Use the `input` parameter to specify its location. It has to be a comma-separated file with at least three columns, and a header row.

- `tag:` A unique tag which represents the target assembly throughout the pipeline. The `tag` and `fasta` file name should not be same, such as `tag.fasta`. This can create file name collisions in the pipeline or result in file overwrite. It is also a good-practice to make all the input files read-only.
- `fasta:` FASTA file for the genome
- `is_masked:` Whether the FASTA is masked or not? Use yes/no to indicate the masking. If the assembly is not masked. The pipeline will soft mask it before annotating it.
- `te_lib [Optional]:` If an assembly is not masked and a TE library is available which cna be used to mask the assembly, the path of the TE library FASTA file can be provided here. If this column is absent and the assembly is not masked, the pipeline will first create a TE library so that it can soft mask the assembly.
- `benchmark [Optional]:` A GFF3 file which can be used to benchmark or compare the results of the pipeline against an existing annotation.

### Advanced inputs for manual resume

If the pipeline fails while processing large datasets, it is advisable to backup the repeat-masked genomes and the BRAKER outputs before attempting a [Nextflow resume](https://www.nextflow.io/docs/latest/cache-and-resume.html#caching-and-resuming). If the resume fails, these outputs from the first pipeline run can be used to setup a manual resume. This can be achieved by providing the repeat-masked genomes under the `fasta` column along with `is_masked` column set to `yes`. The BRAKER outputs can be provided under the following columns,

- `braker_gff3 [Optional]:` BRAKER GFF3 file
- `braker_hints [Optional]:` BRAKER hints file in GFF3 format

The pipeline will automatically skip the repeat modelling, masking and BRAKER steps. It will still perform these steps for those genomes for which these files are not provided. These files are not saved by the pipeline by default. To save the files, set the `repeatmasker_save_outputs` and `braker_save_outputs` parameters to `true`.

## Protein evidence

> ✅ Mandatory `--protein_evidence`

Protein evidence can be provided in two ways. First, a single FASTA file. Second, a list of FASTA files listed in a plain text file. The extension of the text file must be `txt`.

### BRAKER workflow

With these two parameters, the pipeline has sufficient inputs to execute the [BRAKER workflow C](https://github.com/Gaius-Augustus/BRAKER/tree/f58479fe5bb13a9e51c3ca09cb9e137cab3b8471?tab=readme-ov-file#overview-of-modes-for-running-braker) (see Figure 4) in which GeneMark-EP+ is trained on protein spliced alignments, then GeneMark-EP+ generates training data for AUGUSTUS which then performs the final gene prediction.

## RNASeq evidence

> ❔ Optional `--rna_evidence`

RNASeq evidence must be provided through a samplesheet in CSV format which has the following columns,

- `sample:` A sample identifier. The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis.
- `file_1:` A SRA ID for paired-end reads or FASTQ or BAM file
- `file_2:` A FASTQ file if `file_1` is also a FASTQ file and provides paired samples.
- `target_assemblies:` A semicolon `;` separated list of assembly tags from the [assemblysheet input](#assemblysheet-input). If `file_1` points to a BAM file, only a single assembly can be listed under `target_assemblies` for that sample. FASTQ data from `file_1` and `file_2` is aligned against each target assembly. BAM data from `file_1` is considered already aligned against the target assembly and is directly fed to BRAKER.

### BRAKER workflow

If RNASeq evidence is provided, the pipeline executes the [BRAKER workflow D](https://github.com/Gaius-Augustus/BRAKER/tree/f58479fe5bb13a9e51c3ca09cb9e137cab3b8471?tab=readme-ov-file#overview-of-modes-for-running-braker) (see Figure 4) in which GeneMark-ETP is trained with both protein and RNASeq evidence and the training data generated by GeneMark-ETP is used to optimise AUGUSTUS for final gene predictions.

### Preprocessing

RNASeq reads provided in FASTQ files are by default trimmed with [FASTP](https://github.com/OpenGene/fastp). No parameters are provided by default. Although, additional parameters can be provided with `--fastp_extra_args` parameter. After trimming, any sample which does not have `10000` reads left is dropped. This threshold can be specified with the `--min_trimmed_reads` parameter. If trimming was already performed ot it is not desirable, it can be skipped by setting the `--fastp_skip` flag to `true`.

Optionally, [SORTMERNA](https://github.com/sortmerna/sortmerna) can be activated by setting the `--remove_ribo_rna` flag to `true`. A default list of rRNA databases is pre-configured and can be seen in the [assets/rrna-db-defaults.txt](../assets/rrna-db-defaults.txt) file. A path to a custom list of databases can be specified by the `--ribo_database_manifest` parameter.

### Alignment

RNASeq evidence provided as FASTQ files is aligned using [STAR](https://github.com/alexdobin/STAR). The default alignment parameters are,

```bash
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand gunzip -c \
--alignIntronMax $star_max_intron_length
```

where `--star_max_intron_length` is a pipeline parameter and its default value is `16000`. In our experience, the performance of BRAKER predictions is fairly sensitive to this parameter and the parameter value should be based on some estimation of the length of introns in the genes of the target _species_. Additional STAR parameters can be specified with `--star_align_extra_args`.

> [!WARNING]
>
> If pre-aligned RNASeq data is provided as a BAM file and the alignment was not performed with `--outSAMstrandField intronMotif` parameter, the pipeline might trough an error.

## Liftoff annotations

> ❔ Optional `--liftoff_annotations`

In addition to gene prediction with BRAKER, the pipeline also enables gene model transfer from one or more reference assemblies to all the target assemblies. The reference assemblies and the associated gene models must be specified through a CSV file with the following two columns,

- `fasta:` Reference assembly genome in a FASTA file
- `gff3:` Reference assembly gene models in a GFF3 file

[LIFTOFF](https://github.com/agshumate/Liftoff) is used for lifting over the models. The default alignment parameters are,

```bash
-exclude_partial \
-copies \
-polish \
-a $liftoff_coverage \
-s $liftoff_identity
```

where `--liftoff_coverage` and `--liftoff_identity` are pipeline parameters and their default value is `0.9`. After the liftoff, the pipeline filters out any model which is marked as `valid_ORF=False` by [LIFTOFF](https://github.com/agshumate/Liftoff). Then, the BRAKER and LIFTOFF annotations are merged together. During this merge, LIFTOFF models are given precedence over BRAKER models. A region where the LIFTOFF model overlaps a BRAKER model, the BRAKER model is dropped.

## EggNOG-mapper DB

> ❔ Optional `--eggnogmapper_db_dir`, `--eggnogmapper_tax_scope`

EggNOG-mapper is used to add functional annotations to the gene models. The EggNOG-mapper database must be downloaded manually before running the pipeline. The database is available at <http://eggnog5.embl.de/#/app/downloads>. The path to the database folder must be provided with the `--eggnogmapper_db_dir` parameter. The pipeline assumes following directory structure for the database path.

```bash
/path/to/db
├── eggnog.db
├── eggnog.taxa.db
├── eggnog.taxa.db.traverse.pkl
├── eggnog_proteins.dmnd
├── mmseqs
│   ├── mmseqs.db
│   ├── mmseqs.db.dbtype
│   ├── mmseqs.db.index
│   ├── mmseqs.db.lookup
│   ├── mmseqs.db.source
│   ├── mmseqs.db_h
│   ├── mmseqs.db_h.dbtype
│   └── mmseqs.db_h.index
└── pfam
    ├── Pfam-A.clans.tsv.gz
    ├── Pfam-A.hmm
    ├── Pfam-A.hmm.h3f
    ├── Pfam-A.hmm.h3i
    ├── Pfam-A.hmm.h3m
    ├── Pfam-A.hmm.h3m.ssi
    ├── Pfam-A.hmm.h3p
    └── Pfam-A.hmm.idmap
```

An appropriate taxonomic scope for the mapper can be specified with `--eggnogmapper_tax_scope` parameter, otherwise, the pipeline uses teh default value of `1` for the taxonomic scope. Common taxonomic scopes are Eukaryota: `2759`, Viridiplantae: `33090`, Archaea: `2157`, Bacteria: `2` and root: `1`. For a comprehensive list of available scopes, see <http://eggnog5.embl.de/#/app/downloads>.

## Orthology inference input

> ❔ Optional `--orthofinder_annotations`

If there are more than one target assemblies, an orthology inference is performed with [ORTHOFINDER](https://github.com/davidemms/OrthoFinder). Additional annotations can be directly provided for the orthology inference with the `--orthofinder_annotations` parameter. This should be the path to a CSV file with following two columns,

- `tag:` A unique tag which represents the annotation. The `tag` and `fasta` file name should not be same, such as `tag.fasta`. This can create file name collisions in the pipeline or result in file overwrite. It is also a good-practice to make all the input files read-only.
- `fasta:` FASTA file containing protein sequences.

## Iso-forms and full intron support

By default the pipeline allows multiple isoforms from BRAKER. This behavior can be changed by setting the `--allow_isoforms` flag to `false`. Moreover, every intron from every model from BRAKER and LIFTOFF must have support from protein or RNASeq evidence. This is enforced with [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA). This requirement can be removed by setting the `--enforce_full_intron_support` flag to `false`. Or, selectively only applying this criterion to BRAKER models by setting the `--filter_liftoff_by_hints` flag to `false`.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run plant-food-research-open/genepal \
  -revision <version> \
  -profile <docker/singularity/.../institute> \
  --input assemblysheet.csv \
  --protein_evidence proteins.faa \
  --outdir <OUTDIR>
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run plant-food-research-open/genepal -revision main -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './assemblysheet.csv'
outdir: './results/'
protein_evidence: './proteins.faa'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull plant-food-research-open/genepal
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [plant-food-research-open/genepal releases page](https://github.com/plant-food-research-open/genepal/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
