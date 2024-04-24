# PANGENE

[![Lint and -stub on Linux/Docker](https://github.com/PlantandFoodResearch/pangene/actions/workflows/test.yml/badge.svg)](https://github.com/PlantandFoodResearch/pangene/actions/workflows/test.yml)

A NextFlow pipeline for pan-genome annotation. It can also be used for annotation of a single genome.

## Pipeline Flowchart

```mermaid
flowchart TD
    subgraph PrepareAssembly [ ]
    TARGET_ASSEMBLIES
    TE_LIBRARIES
    FASTA_VALIDATE
    fasta_file_from_fasta_validate
    EDTA
    REPEATMODELER
    te_lib_absent_node
    REPEATMASKER
    end

    TARGET_ASSEMBLIES(["[target_assemblies]"])
    TE_LIBRARIES(["[te_libs]"])
    TARGET_ASSEMBLIES --> FASTA_VALIDATE
    FASTA_VALIDATE --- |Fasta|fasta_file_from_fasta_validate(( ))
    fasta_file_from_fasta_validate --> |or|EDTA
    fasta_file_from_fasta_validate --> |default|REPEATMODELER
    REPEATMODELER --- te_lib_absent_node(( ))
    EDTA --- te_lib_absent_node
    TE_LIBRARIES --> REPEATMASKER
    te_lib_absent_node --> REPEATMASKER

    subgraph Fastq [ ]
    FASTQ
    CAT_FASTQ
    FASTQC
    FASTP
    FASTP_FASTQC
    SORTMERNA
    fasta_file_for_star
    STAR
    SAMTOOLS_CAT
    end

    FASTQ([fastq])
    FASTQ --> |Tech. reps|CAT_FASTQ
    CAT_FASTQ --> FASTQC
    FASTQ --> FASTQC
    FASTQC --> FASTP
    FASTP --> FASTP_FASTQC[FASTQC]
    FASTP_FASTQC --> SORTMERNA
    fasta_file_for_star(( ))
    fasta_file_for_star --> |Fasta|STAR
    SORTMERNA --> STAR
    STAR --> SAMTOOLS_CAT

    subgraph Annotation [ ]
    anno_fasta(( ))
    anno_masked_fasta(( ))
    anno_bam(( ))
    EXTERNAL_PROTEIN_SEQS(["[ext_prots]"])
    XREF_ANNOTATIONS(["[xref_annotations]"])
    CAT
    BRAKER3
    GFFREAD
    LIFTOFF
    end

    PrepareAssembly --> |Fasta, Masked fasta|Annotation
    Fastq --> |RNASeq bam|Annotation

    XREF_ANNOTATIONS --> |xref_gff|GFFREAD
    XREF_ANNOTATIONS --> |xref_fasta|LIFTOFF
    GFFREAD --> LIFTOFF
    anno_fasta --> |Fasta|LIFTOFF

    EXTERNAL_PROTEIN_SEQS --> CAT
    anno_masked_fasta --> |Masked fasta|BRAKER3
    anno_bam --> |RNASeq bam|BRAKER3
    CAT --> BRAKER3

    style Fastq fill:#00FFFF21,stroke:#00FFFF21
    style PrepareAssembly fill:#00FFFF21,stroke:#00FFFF21
    style Annotation fill:#00FFFF21,stroke:#00FFFF21
```

## Alpha Release

This release is not fully documented and under alpha testing by the Bioinformatics Team. There are several [outstanding issues](https://github.com/PlantandFoodResearch/pangene/issues) which will be addressed before a general release.

## Plant&Food Users

Download the pipeline to your `/workspace/$USER` folder. Change the parameters defined in the [pfr/params.json](./pfr/params.json) file. Submit the pipeline to SLURM for execution. For a description of the parameters, see [parameters](./docs/parameters.md).

```bash
sbatch ./pfr_pangene
```

## Credits

plantandfoodresearch/pangene scripts were originally written by Jason Shiller. Usman Rashid wrote the NextFLow pipeline.

We thank the following people for their extensive assistance in the development of this pipeline.

- Cecilia Deng [@CeciliaDeng](https://github.com/CeciliaDeng)
- Charles David [@charlesdavid](https://github.com/charlesdavid)
- Chen Wu [@christinawu2008](https://github.com/christinawu2008)
- Leonardo Salgado [@leorippel](https://github.com/leorippel)
- Ross Crowhurst [@rosscrowhurst](https://github.com/rosscrowhurst)
- Susan Thomson [@cflsjt](https://github.com/cflsjt)
- Ting-Hsuan Chen [@ting-hsuan-chen](https://github.com/ting-hsuan-chen)

The pipeline uses nf-core modules contributed by following authors.

<a href="https://github.com/drpatelh"><img src="https://github.com/drpatelh.png" width="50" height="50"></a>
<a href="https://github.com/edmundmiller"><img src="https://github.com/edmundmiller.png" width="50" height="50"></a>
<a href="https://github.com/erikrikarddaniel"><img src="https://github.com/erikrikarddaniel.png" width="50" height="50"></a>
<a href="https://github.com/ewels"><img src="https://github.com/ewels.png" width="50" height="50"></a>
<a href="https://github.com/felixkrueger"><img src="https://github.com/felixkrueger.png" width="50" height="50"></a>
<a href="https://github.com/friederikehanssen"><img src="https://github.com/friederikehanssen.png" width="50" height="50"></a>
<a href="https://github.com/gallvp"><img src="https://github.com/gallvp.png" width="50" height="50"></a>
<a href="https://github.com/grst"><img src="https://github.com/grst.png" width="50" height="50"></a>
<a href="https://github.com/jemten"><img src="https://github.com/jemten.png" width="50" height="50"></a>
<a href="https://github.com/jfy133"><img src="https://github.com/jfy133.png" width="50" height="50"></a>
<a href="https://github.com/joseespinosa"><img src="https://github.com/joseespinosa.png" width="50" height="50"></a>
<a href="https://github.com/kevinmenden"><img src="https://github.com/kevinmenden.png" width="50" height="50"></a>
<a href="https://github.com/kherronism"><img src="https://github.com/kherronism.png" width="50" height="50"></a>
<a href="https://github.com/mashehu"><img src="https://github.com/mashehu.png" width="50" height="50"></a>
<a href="https://github.com/matthdsm"><img src="https://github.com/matthdsm.png" width="50" height="50"></a>
<a href="https://github.com/praveenraj2018"><img src="https://github.com/praveenraj2018.png" width="50" height="50"></a>
<a href="https://github.com/robsyme"><img src="https://github.com/robsyme.png" width="50" height="50"></a>
<a href="https://github.com/toniher"><img src="https://github.com/toniher.png" width="50" height="50"></a>
<a href="https://github.com/vagkaratzas"><img src="https://github.com/vagkaratzas.png" width="50" height="50"></a>

## Citations

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
