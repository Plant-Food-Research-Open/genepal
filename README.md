# PANGENE

A NextFlow pipeline for pan-genome annotation.

## Pipeline Flowchart

```mermaid
flowchart TD
    subgraph PrepareAssembly [ ]
    TARGET_ASSEMBLIES
    TE_LIBRARIES
    FASTA_VALIDATE
    fasta_file_from_fasta_validate
    EDTA
    REPEATMASKER
    end

    TARGET_ASSEMBLIES(["[target_assemblies]"])
    TE_LIBRARIES(["[te_libs]"])
    TARGET_ASSEMBLIES --> FASTA_VALIDATE
    FASTA_VALIDATE --> |Fasta|fasta_file_from_fasta_validate(( ))
    fasta_file_from_fasta_validate --> EDTA
    TE_LIBRARIES --> REPEATMASKER
    EDTA --> |te_lib absent|REPEATMASKER

    subgraph Samplesheet [ ]
    SAMPLESHEET
    CAT_FASTQ
    FASTQC
    FASTP
    FASTP_FASTQC
    SORTMERNA
    fasta_file_for_star
    STAR
    SAMTOOLS_CAT
    end

    SAMPLESHEET([samplesheet])
    SAMPLESHEET --> |Tech. reps|CAT_FASTQ
    CAT_FASTQ --> FASTQC
    SAMPLESHEET --> FASTQC
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
    Samplesheet --> |RNASeq bam|Annotation

    XREF_ANNOTATIONS --> |xref_gff|GFFREAD
    XREF_ANNOTATIONS --> |xref_fasta|LIFTOFF
    GFFREAD --> LIFTOFF
    anno_fasta --> |Fasta|LIFTOFF

    EXTERNAL_PROTEIN_SEQS --> CAT
    anno_masked_fasta --> |Masked fasta|BRAKER3
    anno_bam --> |RNASeq bam|BRAKER3
    CAT --> BRAKER3

    style Samplesheet fill:#00FFFF21,stroke:#00FFFF21
    style PrepareAssembly fill:#00FFFF21,stroke:#00FFFF21
    style Annotation fill:#00FFFF21,stroke:#00FFFF21
```

## Plant&Food Users

Configure the pipeline by modifying `nextflow.config` and submit to SLURM for execution.

```bash
sbatch ./pan_gene_pfr.sh
```

## Third-party Sources

Some software components of this pipeline have been adopted from following third-party sources:

1. nf-core [MIT](https://github.com/nf-core/modules/blob/master/LICENSE): https://github.com/nf-core/modules

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

2. nf-core/rnaseq [MIT](https://github.com/nf-core/rnaseq/blob/master/LICENSE): https://github.com/nf-core/rnaseq
3. rewarewaannotation [MIT](https://github.com/kherronism/rewarewaannotation/blob/master/LICENSE): https://github.com/kherronism/rewarewaannotation
4. assembly_qc [GPL-3.0](https://github.com/Plant-Food-Research-Open/assembly_qc/blob/main/LICENSE): https://github.com/Plant-Food-Research-Open/assembly_qc
