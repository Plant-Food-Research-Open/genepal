---
title: "Braker3_benchmarking"
author: "Jason Shiller"
format: html
editor: visual
---

## Testing EASEL gene prediction and annotation pipeline

# Prep files and directories.

```{r}
Sys.setenv(hdir = "/workspace/hrtjbs/Annotation_Pipeline_2023",
           gdir = "/workspace/hrtjbs/Annotation_Pipeline_2023/genome",
           rnaseqdir = "/input/genomic/plant/Actinidia/chinensis/hort16A/Transcriptome/Illumina/PairEnd/1505KHS-0090",
           oref = "/workspace/cflthc/scratch/2022_Actinidia_TE/08.02_RepeatMasker/CK/CK6901M.chromosomes.only.fa.masked",
           tinyref = "/workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chr6.masked.fa",
           reads = "/workspace/hrtjbs/Annotation_Pipeline_2023/rnaseq",
           trim_reads = "/workspace/hrtjbs/Annotation_Pipeline_2023/rnaseq/trim",
           trinitydir = "/workspace/hrtjbs/Annotation_Pipeline_2023/transcriptAssembly",
           pasa_dir = "/workspace/hrtjbs/Annotation_Pipeline_2023/pasa",
           braker_dir = "/workspace/hrtjbs/Annotation_Pipeline_2023/braker",
           mapping_dir = "/workspace/hrtjbs/Annotation_Pipeline_2023/mapping",
           oredref = "/workspace/cflthc/scratch/2022_Actinidia_TE/08.02_RepeatMasker/RE/Red5.chromosomes.only.fa.masked",
           redref = "/workspace/hrtjbs/Annotation_Pipeline_2023/genome/Red5.chromosomes.only.masked.fa")
```

##Bench test

# First we are using Chromomsome 6 of the CK6901M genome and a subset of RNA-seq samples including shoot, root,leaf,fruit,flower and cane

# Extract Chr6 of CK6901M

```{bash}

module load samtools
mkdir -p $hdir
mkdir -p $gdir
mkdir -p $reads
mkdir -p $trim_reads
#copy genome to working dir and renome slightly
cp $oref /workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chromosomes.only.masked.fa
samtools faidx /workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chromosomes.only.masked.fa
samtools faidx /workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chromosomes.only.masked.fa  chr6 > $tinyref

```

#Get RNA-seq reads

```{bash}
for i in `ls $rnaseqdir | grep -E 'cane3|Flower11|Root1_|shoot5|Fruit_T115|Sink10'`
do cp  $rnaseqdir/$i $reads/
done
```

#Trim reads

```{bash}

cd $hdir
readarray -t r1s < <(ls $reads/*1.fastq.gz)
readarray -t r2s < <(ls $reads/*2.fastq.gz)
readarray -t sn < <(ls  $reads/*1.fastq.gz | awk -F "/" '{print $NF}' | awk -F "_" '{print$2}')

for ((i=0;i<${#r1s[@]};i++));
do

cat << EOF >  ${sn[$i]}.trim.sl
#!/bin/bash -e
#SBATCH -J ${sn[$i]}
#SBATCH --output=${sn[$i]}.trim.out
#SBATCH --error=${sn[$i]}.trim.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00

module load fastp/0.23.2

fastp -l 20 -w 8 -i ${r1s[$i]} \
-I  ${r2s[$i]} -o $trim_reads/${sn[$i]}.fastp.R1.fq.gz -O $trim_reads/${sn[$i]}.fastp.R2.fq.gz
EOF
sbatch ${sn[$i]}.trim.sl
done

```

#Generate de-novo transcriptome assembly with Trinity

```{bash}

mkdir -p $trinitydir
cd $trinitydir


readarray -t r1trimmed < <(ls /workspace/hrtjbs/Annotation_Pipeline_2023/rnaseq/trim/*.R1.fq.gz )
readarray -t r2trimmed < <(ls /workspace/hrtjbs/Annotation_Pipeline_2023/rnaseq/trim/*.R2.fq.gz )

cat << EOF > trinityDenovo.sl
#!/bin/bash -e
#SBATCH -J TDN
#SBATCH --output=TrinityDn.out
#SBATCH --error=TrinityDn.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=7-00:00:00

module load  trinityrnaseq/2.14.0

trinityrnaseq.v2.14.0.sif \
--seqType fq \
--left  `echo ${r1trimmed[@]} | tr ' ' ',' | sed 's/,$//g'` \
--right `echo ${r2trimmed[@]} | tr ' ' ',' | sed 's/,$//g'` \
--max_memory 96G \
--SS_lib_type RF \
--CPU 16 \
--min_contig_length 250 \
--output $trinitydir/trinity.out
EOF
sbatch trinityDenovo.sl


```

#Clean trinity transcripts with Evigene

```{bash}
cd $hdir

mkdir -p clean_trans
cd clean_trans

cat << 'EOF' > Evigene.sl
#!/bin/bash -e
#SBATCH -J eviG
#SBATCH --output=eviG.out
#SBATCH --error=eviG.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=08:00:00

module load  evigene/30-07-2021
tr2aacds.pl -tidy -NCPU=24 -MAXMEM=24000 -log -debug  -cdna /workspace/hrtjbs/Annotation_Pipeline_2023/transcriptAssembly/trinity.out.Trinity.fasta
EOF
sbatch Evigene.sl
```

#prepare congfigs and directory structure for pasa

```{bash}

runID="R1_chr6"

mkdir -p $pasa_dir/$runID\_tmp
mkdir -p $pasa_dir/$runID\_work

cd $pasa_dir/$runID\_work

#copy the input files required

cp /workspace/hrtjbs/Annotation_Pipeline_2023/clean_trans/okayset/trinity.out.Trinity.okay.out.Trinity.cds clean.transcripts.fasta
cp $tinyref genome.fasta

# make the sqlite3 file
touch $runID.db

#edit the conf files to point to the sql db

sed  "s|DATABASE=.*|DATABASE=$pasadir/$runID\_work/$runID.db|g" /workspace/hrtjbs/github_pasa/work_Neo_R2/conf.txt > conf.txt
sed  "s|DATABASE=.*|DATABASE=$pasadir/$runID\_work/$runID.db|g" /workspace/hrtjbs/github_pasa/work_Neo_R2/pasa.annotationCompare.conf.txt > pasa.annotationCompare.conf.txt
sed  "s|DATABASE=.*|DATABASE=$pasadir/$runID\_work/$runID.db|g" /workspace/hrtjbs/github_pasa/work_Neo_R2/pasa.alignAssembly.conf.txt > pasa.alignAssembly.conf.txt



cat conf.txt
cat pasa.annotationCompare.conf.txt
cat pasa.alignAssembly.conf.txt


```

#Run Pasa aligner (to map clean trinity transcripts to genome)

```{bash}
cd /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work
cat << EOF > PASA.sl
#!/bin/bash -e
#SBATCH -J pasa
#SBATCH --output=PASAalign.out
#SBATCH --error=PASAalign.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=23:00:00

#align the seqs
singularity exec -B /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_tmp:/workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_tmp \
-B /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work:/workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work \
/workspace/hrtjbs/github_pasa/pasapipeline_latest.sif \
bash -c 'cd /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work && /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
--MAX_INTRON_LENGTH 25000 \
-c conf.txt \
--TRANSDECODER \
-C \
-R \
--ALIGNER gmap,blat \
--CPU 24 \
-g genome.fasta \
--stringent_alignment_overlap 30.0 \
-t clean.transcripts.fasta'
EOF

sbatch PASA.sl
```

#Genrate training prots to be used by Braker3

```{bash}

cd /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work


singularity exec -B /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_tmp:/workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_tmp \
-B /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work:/workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work \
/workspace/hrtjbs/github_pasa/pasapipeline_latest.sif \
bash -c 'cd /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work  && /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work/R1_chr6.db.assemblies.fasta --pasa_transcripts_gff3 /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work/R1_chr6.db.pasa_assemblies.gff3'

```

#Get the complete proteins only

```{bash}

mkdir -p /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots
cd  /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots
module load samtools
samtools faidx /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work/R1_chr6.db.assemblies.fasta.transdecoder.pep
grep "complete" /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work/R1_chr6.db.assemblies.fasta.transdecoder.pep | awk '{print$1}' | sed s'/>//g'  > /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots/keep.me.pasa
xargs samtools faidx /workspace/hrtjbs/Annotation_Pipeline_2023/pasa/R1_chr6_work/R1_chr6.db.assemblies.fasta.transdecoder.pep < /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots/keep.me.pasa > /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots/pasa.complete.pep.faa

```

#Mapo RNA-seq reads too re genome (first generate star index)

```{bash}
mkdir -p $mapping_dir
cd $mapping_dir
cat << EOF >  StarGenerate.sl
#!/bin/bash -e
#SBATCH -J generate
#SBATCH --output=generate.out
#SBATCH --error=generate.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00

module load  STAR/2.7.10a
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeSAindexNbases 11 \
--genomeDir /workspace/hrtjbs/Annotation_Pipeline_2023/genome/ \
--genomeFastaFiles /workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chr6.masked.fa
EOF
sbatch StarGenerate.sl
```

# Map RNA-seq reads to genome (note outSAMstrandField 'intronMotif' is crucial for braker3)

```{bash}
mkdir -p $mapping_dir
cd $mapping_dir
readarray -t r1trimmed < <(ls /workspace/hrtjbs/Annotation_Pipeline_2023/rnaseq/trim/*.R1.fq.gz )
readarray -t r2trimmed < <(ls /workspace/hrtjbs/Annotation_Pipeline_2023/rnaseq/trim/*.R2.fq.gz )
readarray -t sn < <(ls /workspace/hrtjbs/Annotation_Pipeline_2023/rnaseq/trim/*.R1.fq.gz | awk -F "." '{print$1}' | awk -F "/" '{print$NF}')

for ((i=0;i<${#r1trimmed[@]};i++));
do

cat << EOF >  ${sn[$i]}.star.sl
#!/bin/bash -e
#SBATCH -J ${sn[$i]}
#SBATCH --output=${sn[$i]}.trim.out
#SBATCH --error=${sn[$i]}.trim.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00

module load  STAR/2.7.10a
cd $mapping_dir
STAR --runThreadN 8 \
--readFilesCommand zcat \
--readFilesIn  ${r1trimmed[$i]}  ${r2trimmed[$i]} \
--genomeDir /workspace/hrtjbs/Annotation_Pipeline_2023/genome \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${sn[$i]} \
--outSAMunmapped Within
EOF
done
```

#Fetch orthodb proteins for use with the braker annotation and merge with the proteins created with trintiy / pasa

```{bash}
cd /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots

wget https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Viridiplantae.fa.gz .
gunzip Viridiplantae.fa.gz
cat Viridiplantae.fa pasa.complete.pep.faa > viridiplanta.odb11.plus.pasa.fa
```

#Run braker 3 (This run includes proteins from orthodb + prots generted with trinity pasa + mapped rna-seq reads )

```{bash}
mkdir -p $braker_dir/R1

cd $braker_dir/R1
cp /workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chr6.masked.fa genome.fasta

#copy over the training prots

cp /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots/viridiplanta.odb11.plus.pasa.fa .

cat << 'EOF' > Braker3_R1.sl
#!/bin/bash -e
#SBATCH -J Braker3
#SBATCH --output=braker3_R1.out
#SBATCH --error=braker3_R1.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01-00:00:00

module load conda
conda deactivate
module load pfr-python3/3.9.13
module load samtools/1.16
module load singularity/3

singularity exec -B /workspace/hrtjbs/Annotation_Pipeline_2023/braker/R1:/workspace/hrtjbs/Annotation_Pipeline_2023/braker/R1 \
-B /workspace/hrtjbs/Annotation_Pipeline_2023/mapping:/workspace/hrtjbs/Annotation_Pipeline_2023/mapping \
-B /workspace/hrtjbs/augustus/:/workspace/hrtjbs/augustus/ \
/workspace/hrtjbs/software/braker3.sif braker.pl \
--AUGUSTUS_CONFIG_PATH=/workspace/hrtjbs/augustus/config \
--species=CK6901M_R1 \
--gff3 \
--genome=genome.fasta \
--prot_seq=viridiplanta.odb11.plus.pasa.fa \
--alternatives-from-evidence=true \
--rnaseq_sets_ids=cane3Aligned.sortedByCoord.out,Flower11Aligned.sortedByCoord.out,FruitAligned.sortedByCoord.out,LeafAligned.sortedByCoord.out,Root1Aligned.sortedByCoord.out,shoot5Aligned.sortedByCoord.out \
--rnaseq_sets_dirs=/workspace/hrtjbs/Annotation_Pipeline_2023/mapping \
--workingdir=/workspace/hrtjbs/Annotation_Pipeline_2023/braker/R1 \
--GENEMARK_PATH=${ETP}/gmes \
--threads 16
EOF

sbatch Braker3_R1.sl
```

#run braker3 again using just the proteins generated with pasa + trinity + mapped rna-seq reads but without no orthodb prots

```{bash}

mkdir -p $braker_dir/R2

cd $braker_dir/R2
cp /workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chr6.masked.fa genome.fasta

#copy over the training prots

cp /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots/pasa.complete.pep.faa .

cat << 'EOF' > Braker3_R2.sl
#!/bin/bash -e
#SBATCH -J Braker3
#SBATCH --output=braker3_R2.out
#SBATCH --error=braker3_R2.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01-00:00:00

module load conda
conda deactivate
module load pfr-python3/3.9.13
module load samtools/1.16
module load singularity/3

singularity exec -B /workspace/hrtjbs/Annotation_Pipeline_2023/braker/R2:/workspace/hrtjbs/Annotation_Pipeline_2023/braker/R2 \
-B /workspace/hrtjbs/Annotation_Pipeline_2023/mapping:/workspace/hrtjbs/Annotation_Pipeline_2023/mapping \
-B /workspace/hrtjbs/augustus/:/workspace/hrtjbs/augustus/ \
/workspace/hrtjbs/software/braker3.sif braker.pl \
--AUGUSTUS_CONFIG_PATH=/workspace/hrtjbs/augustus/config \
--species=CK6901M_R2 \
--gff3 \
--genome=genome.fasta \
--prot_seq=pasa.complete.pep.faa \
--alternatives-from-evidence=true \
--rnaseq_sets_ids=cane3Aligned.sortedByCoord.out,Flower11Aligned.sortedByCoord.out,FruitAligned.sortedByCoord.out,LeafAligned.sortedByCoord.out,Root1Aligned.sortedByCoord.out,shoot5Aligned.sortedByCoord.out \
--rnaseq_sets_dirs=/workspace/hrtjbs/Annotation_Pipeline_2023/mapping \
--workingdir=/workspace/hrtjbs/Annotation_Pipeline_2023/braker/R2 \
--GENEMARK_PATH=${ETP}/gmes \
--threads 16
EOF

sbatch Braker3_R2.sl
```

#Run braker3 again using just the mapped rna-seq reads and the orthodb prots

```{bash}
mkdir -p $braker_dir/R3

cd $braker_dir/R3
cp /workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chr6.masked.fa genome.fasta

copy over the training prots

cp /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots/Viridiplantae.fa .

cat << 'EOF' > Braker3_R3.sl
#!/bin/bash -e
#SBATCH -J Braker3
#SBATCH --output=braker3_R3.out
#SBATCH --error=braker3_R3.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00

module load conda
conda deactivate
module load pfr-python3/3.9.13
module load samtools/1.16
module load singularity/3

singularity exec -B /workspace/hrtjbs/Annotation_Pipeline_2023/braker/R3:/workspace/hrtjbs/Annotation_Pipeline_2023/braker/R3 \
-B /workspace/hrtjbs/Annotation_Pipeline_2023/mapping:/workspace/hrtjbs/Annotation_Pipeline_2023/mapping \
-B /workspace/hrtjbs/augustus/:/workspace/hrtjbs/augustus/ \
/workspace/hrtjbs/software/braker3.sif braker.pl \
--AUGUSTUS_CONFIG_PATH=/workspace/hrtjbs/augustus/config \
--species=CK6901M_R3_2 \
--gff3 \
--genome=genome.fasta \
--prot_seq=Viridiplantae.fa \
--alternatives-from-evidence=true \
--rnaseq_sets_ids=cane3Aligned.sortedByCoord.out,Flower11Aligned.sortedByCoord.out,FruitAligned.sortedByCoord.out,LeafAligned.sortedByCoord.out,Root1Aligned.sortedByCoord.out,shoot5Aligned.sortedByCoord.out \
--rnaseq_sets_dirs=/workspace/hrtjbs/Annotation_Pipeline_2023/mapping \
--workingdir=/workspace/hrtjbs/Annotation_Pipeline_2023/braker/R3 \
--threads 16
EOF

sbatch Braker3_R3.sl
```

#Run braker 3 including kiwifruit prots + orthodb prots + trin pasa prots + mapped rna-seq reads
To test a method which uses our current kiwifruit resources to better annotate new kiwifruit genomes. Created a file with all our kiwifruit proteins and use this in conjunction with orthodb prots

```{bash}

mkdir -p $braker_dir/R4
cd $braker_dir/R4

#copy over training proteins (red5 + russell + orthodb)
cp /workspace/hrtjbs/Annotation_Pipeline_2023/genome/CK6901M.chr6.masked.fa genome.fasta
cat /output/genomic/fairGenomes/Plant/Actinidia/chinensis/var_chinensis/male/2x/assembly_russell/v2.1/RU01.20221117162301.primary.pep.fasta /output/genomic/fairGenomes/Plant/Actinidia/chinensis/var_chinensis/female/2x/assembly_red5/v2/MA20.20220711141824.pep.fasta /workspace/hrtjbs/Annotation_Pipeline_2023/trainingprots/viridiplanta.odb11.plus.pasa.fa > red5_russell_orthodb.fa

cat << 'EOF' > Braker3_R4.sl
#!/bin/bash -e
#SBATCH -J Braker3
#SBATCH --output=braker3_R4.out
#SBATCH --error=braker3_R4.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00

module load conda
conda deactivate
module load pfr-python3/3.9.13
module load samtools/1.16
module load singularity/3

singularity exec -B /workspace/hrtjbs/Annotation_Pipeline_2023/braker/R4:/workspace/hrtjbs/Annotation_Pipeline_2023/braker/R4 \
-B /workspace/hrtjbs/Annotation_Pipeline_2023/mapping:/workspace/hrtjbs/Annotation_Pipeline_2023/mapping \
-B /workspace/hrtjbs/augustus/:/workspace/hrtjbs/augustus/ \
/workspace/hrtjbs/software/braker3.sif braker.pl \
--AUGUSTUS_CONFIG_PATH=/workspace/hrtjbs/augustus/config \
--species=CK6901M_R3_4 \
--gff3 \
--genome=genome.fasta \
--prot_seq=red5_russell_orthodb.fa \
--alternatives-from-evidence=true \
--rnaseq_sets_ids=cane3Aligned.sortedByCoord.out,Flower11Aligned.sortedByCoord.out,FruitAligned.sortedByCoord.out,LeafAligned.sortedByCoord.out,Root1Aligned.sortedByCoord.out,shoot5Aligned.sortedByCoord.out \
--rnaseq_sets_dirs=/workspace/hrtjbs/Annotation_Pipeline_2023/mapping \
--workingdir=/workspace/hrtjbs/Annotation_Pipeline_2023/braker/R4 \
--threads 16
EOF

sbatch Braker3_R4.sl
```
