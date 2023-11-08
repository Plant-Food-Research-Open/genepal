# PAN-GENE
A NextFlow pipeline for pan-genome annotation.

## Pipeline Flowchart

```mermaid
flowchart TD
    ribo_db((ribo_db))
    SAMPLESHEET((samples))
    TE_LIBRARIES(("[te_libs]"))
    TARGET_ASSEMBLIES(("[assemblies]"))
    EXTERNAL_PROTEIN_SEQS(("[ext_prots]"))
    
    GUNZIP_PROT[GUNZIP]
    GUNZIP_TE[GUNZIP]
    SKIP_EDTA{Skip EDTA}
    pend((dev))
    
    TE_LIBRARIES --> GUNZIP_TE
    GUNZIP_TE --> SKIP_EDTA
    
    TARGET_ASSEMBLIES --> GUNZIP
    GUNZIP --> FASTA_VALIDATE
    FASTA_VALIDATE --> FASTA_PERFORM_EDTA
    FASTA_VALIDATE --> SKIP_EDTA
    
    SKIP_EDTA --> REPEATMASKER
    FASTA_PERFORM_EDTA --> REPEATMASKER
    REPEATMASKER --> STAR_GENOMEGENERATE

    SAMPLESHEET --> SAMPLESHEET_CHECK
    SAMPLESHEET_CHECK --> |Technical replicates|CAT_FASTQ
    CAT_FASTQ --> FASTQC
    SAMPLESHEET_CHECK --> FASTQC
    FASTQC --> FASTP
    
    ribo_db --> SORTMERNA
    FASTP --> SORTMERNA
    SORTMERNA --> STAR_ALIGN
    STAR_GENOMEGENERATE --> STAR_ALIGN
    STAR_ALIGN --> GROUP_BY_ASSEMBLY([Group by assembly])
    GROUP_BY_ASSEMBLY --> SAMTOOLS_CAT
    SAMTOOLS_CAT --> |RNASeq bam|BRAKER3

    REPEATMASKER --> BRAKER3

    EXTERNAL_PROTEIN_SEQS --> GUNZIP_PROT
    GUNZIP_PROT --> CAT
    CAT --> BRAKER3
    
    BRAKER3 --> pend

    subgraph Params
    TARGET_ASSEMBLIES
    TE_LIBRARIES
    SAMPLESHEET
    ribo_db
    EXTERNAL_PROTEIN_SEQS
    end

    subgraph GenomePrep
    GUNZIP
    FASTA_VALIDATE
    GUNZIP_TE
    FASTA_PERFORM_EDTA
    SKIP_EDTA
    REPEATMASKER
    STAR_GENOMEGENERATE
    end

    subgraph Braker
    CAT
    GUNZIP_PROT
    BRAKER3
    end

    subgraph SamplePrep
    SAMPLESHEET_CHECK
    CAT_FASTQ
    FASTQC
    FASTP
    SORTMERNA
    STAR_ALIGN
    GROUP_BY_ASSEMBLY
    SAMTOOLS_CAT
    end

    style Params fill:#00FFFF21,stroke:#00FFFF21
    style GenomePrep fill:#00FFFF21,stroke:#00FFFF21
    style SamplePrep fill:#00FFFF21,stroke:#00FFFF21
    style Braker fill:#00FFFF21,stroke:#00FFFF21
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