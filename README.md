# plant-food-research-open/genepal

[![GitHub Actions CI Status](https://github.com/plant-food-research-open/genepal/actions/workflows/ci.yml/badge.svg)](https://github.com/plant-food-research-open/genepal/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/plant-food-research-open/genepal/actions/workflows/linting.yml/badge.svg)](https://github.com/plant-food-research-open/genepal/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.14195006-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.14195006)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda ❌](http://img.shields.io/badge/run%20with-conda%20❌-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/plant-food-research-open/genepal)

## Introduction

**plant-food-research-open/genepal** is a bioinformatics pipeline for single genome, phased genomes and pan-genome annotation. An overview is shown in the [Pipeline Flowchart](#pipeline-flowchart) and the references for the tools are listed in [CITATIONS.md](./CITATIONS.md). Protein coding gene structures are predicted with [BRAKER](https://github.com/Gaius-Augustus/BRAKER) which uses GeneMark-ES/ET/EP+/ETP. These tools require a [license](https://genemark.bme.gatech.edu/license_download.cgi) for commercial works.

## Pipeline Flowchart

<p align="center"><img src="docs/img/genepal.png"></p>

- [fasta_validator](https://github.com/linsalrob/fasta_validator): Validate genome FASTA
- [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) or [EDTA](https://github.com/oushujun/EDTA): Create TE library
- [RepeatMasker](https://github.com/rmhubley/RepeatMasker): Soft mask the genome fasta
- [sra-tools](https://github.com/ncbi/sra-tools): RNASeq data download from SRA
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc), [fastp](https://github.com/OpenGene/fastp), [SortMeRNA](https://github.com/sortmerna/sortmerna): QC, trim and filter RNASeq evidence
- [STAR](https://github.com/alexdobin/STAR): RNASeq alignment
- [cat](https://github.com/coreutils/coreutils/blob/master/src/cat.c): Concatenate protein FASTA files
- [BRAKER](https://github.com/Gaius-Augustus/BRAKER): Predict protein coding gene structures with GeneMark-ES/ET/EP/ETP and AUGUSTUS
  - Directly provided BAM files should be `--outSAMstrandField intronMotif` compliant
  - With protein evidence alone, [BRAKER workflow C](https://github.com/Gaius-Augustus/BRAKER/tree/f58479fe5bb13a9e51c3ca09cb9e137cab3b8471?tab=readme-ov-file#overview-of-modes-for-running-braker) is executed
  - With protein plus RNASeq evidence, [BRAKER workflow D](https://github.com/Gaius-Augustus/BRAKER/tree/f58479fe5bb13a9e51c3ca09cb9e137cab3b8471?tab=readme-ov-file#overview-of-modes-for-running-braker) is executed
- [Liftoff](https://github.com/agshumate/Liftoff): Optionally, liftoff annotations from reference genome FASTA/GFF
- [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA): Optionally, ensure that each BRAKER or both BRAKER and Liftoff models have [full intron support](./docs/usage.md#iso-forms-and-full-intron-support)
- [AGAT](https://github.com/NBISweden/AGAT)
  - Merge multi-reference liftoffs
  - Remove liftoff transcripts marked by _valid_ORF=False_
  - Remove liftoff genes with any intron shorter than 10 bp
  - Remove rRNA and tRNA from liftoff
  - Optionally, allow or remove iso-forms
  - Remove BRAKER models from Liftoff loci
  - Merge Liftoff and BRAKER models
  - Optionally, remove models without any EggNOG-mapper hits
- [EggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper): Add functional annotation to gff
- [GenomeTools](https://github.com/genometools/genometools): GFF format validation
- [GffRead](https://github.com/gpertea/gffread): Extraction of protein sequences
- [OrthoFinder](https://github.com/davidemms/OrthoFinder): Perform phylogenetic orthology inference across genomes
- [GffCompare](https://github.com/gpertea/gffcompare): Compare and benchmark against an existing annotation
- [BUSCO](https://gitlab.com/ezlab/busco): Completeness statistics for genome and annotation through proteins
- [R Markdown](https://rmarkdown.rstudio.com): Specialized pangene analysis
- [MultiQC](https://docs.seqera.io/multiqc): Exhaustive QC statistics

## Usage

Refer to [usage](./docs/usage.md), [parameters](./docs/parameters.md) and [output](./docs/output.md) documents for details.

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare an assemblysheet with your input genomes that looks as follows:

`assemblysheet.csv`:

```csv
tag         ,fasta              ,is_masked
a_thaliana  ,/path/to/genome.fa ,yes
```

Each row represents an input genome and the fields are:

- `tag:` A unique tag which represents the genome throughout the pipeline
- `fasta:` fasta file for the genome
- `is_masked`: yes or no to denote whether the fasta file is already masked or not

#### `--min_contig_length`
- **Description**: Minimum length (in base pairs) of contigs to include in the analysis.
- **Default**: 5000
- **Example**:
    ```bash
    nextflow run main.nf --min_contig_length 10000
    ```
    This will exclude all contigs shorter than 10,000 bp from the analysis.


At minimum, a file with proteins as evidence is also required. Now, you can run the pipeline using:

```bash
nextflow run plant-food-research-open/genepal \
  -revision <version> \
  -profile <docker/singularity/.../institute> \
  --input assemblysheet.csv \
  --protein_evidence proteins.faa \
  --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Plant&Food Users

Download the pipeline to your `/workspace/$USER` folder. Change the parameters defined in the [pfr/params.json](./pfr/params.json) file. Submit the pipeline to SLURM for execution.

```bash
sbatch ./pfr_genepal
```

## Credits

plant-food-research-open/genepal workflows were originally scripted by Jason Shiller ([@jasonshiller](https://github.com/jasonshiller)). Usman Rashid ([@gallvp](https://github.com/gallvp)) wrote the Nextflow pipeline.

We thank the following people for their extensive assistance in the development of this pipeline:

- Cecilia Deng [@CeciliaDeng](https://github.com/CeciliaDeng)
- Charles David [@charlesdavid](https://github.com/charlesdavid)
- Chen Wu [@christinawu2008](https://github.com/christinawu2008)
- Leonardo Salgado [@leorippel](https://github.com/leorippel)
- Ross Crowhurst [@rosscrowhurst](https://github.com/rosscrowhurst)
- Susan Thomson [@cflsjt](https://github.com/cflsjt)
- Ting-Hsuan Chen [@ting-hsuan-chen](https://github.com/ting-hsuan-chen)

The pipeline uses nf-core modules contributed by following authors:

<a href="https://github.com/gallvp"><img src="https://github.com/gallvp.png" width="50" height="50"></a>
<a href="https://github.com/drpatelh"><img src="https://github.com/drpatelh.png" width="50" height="50"></a>
<a href="https://github.com/jfy133"><img src="https://github.com/jfy133.png" width="50" height="50"></a>
<a href="https://github.com/joseespinosa"><img src="https://github.com/joseespinosa.png" width="50" height="50"></a>
<a href="https://github.com/kevinmenden"><img src="https://github.com/kevinmenden.png" width="50" height="50"></a>
<a href="https://github.com/toniher"><img src="https://github.com/toniher.png" width="50" height="50"></a>
<a href="https://github.com/matthdsm"><img src="https://github.com/matthdsm.png" width="50" height="50"></a>
<a href="https://github.com/grst"><img src="https://github.com/grst.png" width="50" height="50"></a>
<a href="https://github.com/friederikehanssen"><img src="https://github.com/friederikehanssen.png" width="50" height="50"></a>
<a href="https://github.com/erikrikarddaniel"><img src="https://github.com/erikrikarddaniel.png" width="50" height="50"></a>
<a href="https://github.com/edmundmiller"><img src="https://github.com/edmundmiller.png" width="50" height="50"></a>
<a href="https://github.com/adamrtalbot"><img src="https://github.com/adamrtalbot.png" width="50" height="50"></a>
<a href="https://github.com/vagkaratzas"><img src="https://github.com/vagkaratzas.png" width="50" height="50"></a>
<a href="https://github.com/kherronism"><img src="https://github.com/kherronism.png" width="50" height="50"></a>
<a href="https://github.com/spficklin"><img src="https://github.com/spficklin.png" width="50" height="50"></a>
<a href="https://github.com/robsyme"><img src="https://github.com/robsyme.png" width="50" height="50"></a>
<a href="https://github.com/priyanka-surana"><img src="https://github.com/priyanka-surana.png" width="50" height="50"></a>
<a href="https://github.com/praveenraj2018"><img src="https://github.com/praveenraj2018.png" width="50" height="50"></a>
<a href="https://github.com/nvnieuwk"><img src="https://github.com/nvnieuwk.png" width="50" height="50"></a>
<a href="https://github.com/muffato"><img src="https://github.com/muffato.png" width="50" height="50"></a>
<a href="https://github.com/maxulysse"><img src="https://github.com/maxulysse.png" width="50" height="50"></a>
<a href="https://github.com/mashehu"><img src="https://github.com/mashehu.png" width="50" height="50"></a>
<a href="https://github.com/mahesh-panchal"><img src="https://github.com/mahesh-panchal.png" width="50" height="50"></a>
<a href="https://github.com/jvhagey"><img src="https://github.com/jvhagey.png" width="50" height="50"></a>
<a href="https://github.com/jemten"><img src="https://github.com/jemten.png" width="50" height="50"></a>
<a href="https://github.com/felixkrueger"><img src="https://github.com/felixkrueger.png" width="50" height="50"></a>
<a href="https://github.com/ewels"><img src="https://github.com/ewels.png" width="50" height="50"></a>
<a href="https://github.com/charles-plessy"><img src="https://github.com/charles-plessy.png" width="50" height="50"></a>
<a href="https://github.com/bunop"><img src="https://github.com/bunop.png" width="50" height="50"></a>
<a href="https://github.com/abhi18av"><img src="https://github.com/abhi18av.png" width="50" height="50"></a>

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use plant-food-research-open/genepal for your analysis, please cite it as:

> **genepal: A Nextflow pipeline for genome and pan-genome annotation.**
>
> Usman Rashid, Jason Shiller, Ross Crowhurst, Chen Wu, Ting-Hsuan Chen, Leonardo Salgado, Charles David, Sarah Bailey, Ignacio Carvajal, Anand Rampadarath, Ken Smith, Liam Le Lievre, Cecilia Deng, Susan Thomson
>
> _zenodo_. 2024. doi: [10.5281/zenodo.14195006](https://doi.org/10.5281/zenodo.14195006).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
