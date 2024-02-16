Pangene_gene.selection
================
Jason Shiller

**Aim:** Take two gene prediction sets for a genome, one generated form
liftoff and the other from Braker3 and output a gtf file that  
contains a new set which contains all the liftoff genes and any
additional Braker3 gene predictions that do not overlap liftoff genes.  
Should also collect information on places where braker3 and Liftoff
agreed

#### Set env variables

``` r
Sys.setenv('hdir' = "/workspace/hrtjbs/pangene_2024/pangene/prototyping", 
           'donorGenome' = "/workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/Russell_V2a.chromosomes.only.fsa",
           'recipientGenome' = '/workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/Red5_V2.chromosomes.only.fsa',
           'donorgff' = '/workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/RU01.20221117162301.no.extras.primary.gff3',
           'brakergff' = '/workspace/hrtjbs/Annotation_Pipeline_2023/braker/RED5_bams_plus_orthoDb_plusPasa/braker.gff3',
           'liftoff' = '/workspace/hrtjbs/software/liftoff/liftoff.sif')
```

#### Liftover genes from Russell to Red5

liftoff.sif = liftoff:1.6.3–pyhdfd78af_0

Remove extra contigs from the Russell assembly

``` bash

grep  -vE '000250F|000252F|000254F|000256F' /workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/RU01.20221117162301.primary.gff3 > /workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/RU01.20221117162301.no.extras.primary.gff3
```

``` bash
mkdir -p $hdir/liftoff
cd $hdir/liftoff

cat << EOF >  liftoff.cp.09.09.sl
#!/bin/bash -e
#SBATCH -J Liftoff
#SBATCH --output=liftoff.%j.out
#SBATCH --error=liftoff.%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00

module load singularity

singularity exec -B /workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff:/workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff \
-B /workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/:/workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/ \
-B /workspace/hrtjbs/Annotation_Pipeline_2023/braker:/workspace/hrtjbs/Annotation_Pipeline_2023/braker \
$liftoff liftoff -g $donorgff \
-p 8 \
-o $hdir/liftoff/red5.from.Russell.cp.09.09.gff3 \
-u $hdir/liftoff/red5.from.Russell.cp.09.09.unmapped.txt \
-exclude_partial \
-copies \
-polish \
-a 0.9 \
-s 0.9 \
$recipientGenome $donorGenome
EOF
sbatch liftoff.cp.09.09.sl
```

Run gffcompare to find overlapping and unique features in each gff3

``` bash
mkdir -p $hdir/gffcompare
cd $hdir/gffcompare

cat << EOF >  gffcompare.cp.09.09.sl
#!/bin/bash -e
#SBATCH -J gffcompare
#SBATCH --output=gffcompare.%j.out
#SBATCH --error=gffcompare.%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00


module load gffcompare/0.12.2
gffcompare -r /workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.gff3_polished -G -o red5.from.Russell.cp.09.09.vs.braker3 $brakergff
EOF
sbatch gffcompare.cp.09.09.sl
```

Examine the results and write to file

``` r
library(dplyr)

classcodes = read.csv("/workspace/hrtjbs/pangene_2024/pangene/prototyping/gffcompare/extended_gffcompare_class_codes.csv")

#funtion to summarise output of gff read
sum.codes = function(orthofinder.output){
  df = read.table(orthofinder.output,sep = "\t")
  df = df %>% group_by(V4) %>%
    summarise(Count = n()) %>% 
  rename( V4 = "Class.Code"  ) %>%
  left_join(y = classcodes,by = 'Class.Code') %>% dplyr::select(Description,Count,Class.Code)
  return(df)
}

#summarise results
red5.from.Russell.cp.09.09.summary.codes = sum.codes("/workspace/hrtjbs/pangene_2024/pangene/prototyping/gffcompare/red5.from.Russell.cp.09.09.vs.braker3.tracking")

write.csv(file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/gffcompare/red5.from.Russell.cp.09.09.summary.codes.csv",x = red5.from.Russell.cp.09.09.summary.codes)
```

\*Note, Codes P ( Single exon genes within 2kb of ref gene) and U
(Intergenic - no overlap) should be retained from the Braker set

``` bash
mkdir -p $hdir/merge_purge

cd $hdir/merge_purge

cat << 'EOF' >  agat.extract.merge.sl
#!/bin/bash -e
#SBATCH -J agat.extract.merge
#SBATCH --output=agat.extract.merge.%j.out
#SBATCH --error=agat.extract.merge.%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00


#extract the IDs of genes which do not overlap liftoff genes (codes u and p )
module load conda
conda activate /workspace/hrtjbs/agat

#awk -F '\t' '$4 == "p" || $4 == "u" { split($5, a, /[|:]/); print a[2] }' /workspace/hrtjbs/pangene_2024/pangene/prototyping/gffcompare/red5.from.Russell.cp.09.09.vs.braker3.tracking | sort | uniq > red5.from.Russell.cp.09.09.vs.braker3.keep

#agat_sp_filter_feature_from_keep_list.pl --gff /workspace/hrtjbs/Annotation_Pipeline_2023/braker/RED5_bams_plus_orthoDb_plusPasa/braker.gff3 \
# --keep_list red5.from.Russell.cp.09.09.vs.braker3.keep \
# --output red5.braker.up.removed.gff3

agat_sp_merge_annotations.pl --gff /workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.gff3 \
--gff red5.braker.up.removed.gff3 \
--out outFile 
EOF

sbatch agat.extract.merge.sl
```

\*Note there were 6524 genes on the list of braker3 genes to be kept
containing a total of 7387 transcripts

Rename outfile (forgot to specify this properly in the script)

``` bash

mv $hdir/merge_purge/outFile  $hdir/merge_purge/braker.liftoff.agat.merge.gff3 
```

Check the output and see gene number etc

``` r
install.packages("GenomicFeatures")
library(GenomicFeatures)
library(rtracklayer)

          
summarise.gff = function(gff_file){
  gff_data = import.gff(gff_file)
  df = as.data.frame(table(mcols(gff_data)$type))
  colnames(df) = c("Feature","Count")
return(df)}
  

  
merged.summary = summarise.gff(gff_file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/braker.liftoff.agat.merge.gff3")
donor.summary = summarise.gff(gff_file = "/workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/RU01.20221117162301.no.extras.primary.gff3")
braker.summary = summarise.gff(gff_file = "/workspace/hrtjbs/Annotation_Pipeline_2023/braker/RED5_bams_plus_orthoDb_plusPasa/braker.gff3")
braker.retained.summary = summarise.gff(gff_file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/red5.braker.up.removed.gff3")
liftoff.summary = summarise.gff(gff_file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.gff3_polished")


write.csv(x = merged.summary,file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/merged.gff.summary.csv")
write.csv(x = donor.summary,file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/donor.gff.summary.csv")
write.csv(x = braker.summary,file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/braker.gff.summary.csv")
write.csv(x = braker.retained.summary,file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/braker.retained.gff.summary.csv")
```

merged summary = 37127 genes (38076 transcripts)  
donor summary = 33826 genes (33826 transcripts)  
braker retained summary = 6524 genes (7387) transcipts

**The expected content of merged file should be 33826 genes + 6524 =
40,350**

Lets try and understand this number:

Firstly: many genes did (3160) not map at all during liftoff, these can
be found in here:  
/workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.unmapped.txt  

**40,350 - 3160 = 37190 (Still doesn’t match)**

The problem seems to be that if one of the isoforms of a gene doesn’t
overlap with a liftover gene but one does the  
gene will be kept in the braker list. This is a problem because during
the merge stage the good isoform and the liftover gene get merged into
the  
name gene even though they don’t overlap.

The way to over come this is to select at the transcript level and then
use agat_convert_sp_gxf2gxf.pl to re format the gff3. This way the
braker gene model gets re-sized  
to just match the size of the good transcript (rather than the original
‘geme’ which overlapped the liftoff gene.

``` bash
cd $hdir/merge_purge

cat << 'EOF' >  agat.extract.merge.v2.sl
#!/bin/bash -e
#SBATCH -J agat.extract.merge
#SBATCH --output=agat.extract.merge.%j.out
#SBATCH --error=agat.extract.merge.%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00


#extract the IDs of genes which do not overlap liftoff genes (codes u and p )

module load conda
conda activate /workspace/hrtjbs/agat

#Note we are extracting the transcript names not gene names

awk -F '\t' '$4 == "p" || $4 == "u" { split($5, a, /[|:]/); print a[3] }' /workspace/hrtjbs/pangene_2024/pangene/prototyping/gffcompare/red5.from.Russell.cp.09.09.vs.braker3.tracking | sort | uniq > red5.from.Russell.cp.09.09.vs.braker3.keep.v2

agat_sp_filter_feature_from_keep_list.pl --gff /workspace/hrtjbs/Annotation_Pipeline_2023/braker/RED5_bams_plus_orthoDb_plusPasa/braker.gff3 \
--keep_list red5.from.Russell.cp.09.09.vs.braker3.keep.v2 \
--output red5.braker.up.removed.v2.gff3

agat_convert_sp_gxf2gxf.pl -g red5.braker.up.removed.v2.gff3 -o red5.braker.up.removed.v2.clean.gff3

agat_sp_merge_annotations.pl --gff /workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.gff3 \
--gff red5.braker.up.removed.v2.clean.gff3 \
--out braker.liftoff.agat.merge.v2.gff3
EOF

sbatch agat.extract.merge.v2.sl
```

Examine the results

``` r
braker.liftoff.agat.merge.v2.gff3 = read.csv("/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/braker.liftoff.agat.merge.v2.gff3",sep = "\t",header = FALSE,skip = 3)

red5.braker.up.removed.v2.clean.gff3 = read.csv("/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/red5.braker.up.removed.v2.clean.gff3",sep = "\t",header = FALSE,skip = 3)
red5.braker.up.removed.v2.clean.gff3  %>% filter(V3 == "mRNA") %>% group_by(V2 )  %>% summarise(Count = n())
```

The numbers still don’t add up. I found the reason. When we select
transcripts it is selecting the gene associated with that transcript as
well as all the other alternative  
transcripts. So the way to solve it is to exclude transcripts form the
file and using;

    agat_sp_filter_feature_from_kill_list.pl
    agat_sp_filter_feature_from_kill_list.pl --gff infile.gff --kill_list file.txt  [ --output outfile ]

We will need to make this list in a slighty different way

``` bash
cd $hdir/merge_purge

cat << 'EOF' >  agat.extract.merge.v3.sl
#!/bin/bash -e
#SBATCH -J agat.extract.merge
#SBATCH --output=agat.extract.merge.%j.out
#SBATCH --error=agat.extract.merge.%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00


#extract the IDs of genes which do not overlap liftoff genes (codes u and p )

module load conda
conda activate /workspace/hrtjbs/agat

#Note we are extracting the transcript names not gene names and now we want the transcripts to exclude not include

awk -F '\t' '$4 != "p" && $4 != "u" { split($5, a, /[|:]/); print a[3] }' /workspace/hrtjbs/pangene_2024/pangene/prototyping/gffcompare/red5.from.Russell.cp.09.09.vs.braker3.tracking | sort | uniq > red5.from.Russell.cp.09.09.vs.braker3.kill

agat_sp_filter_feature_from_kill_list.pl --gff /workspace/hrtjbs/Annotation_Pipeline_2023/braker/RED5_bams_plus_orthoDb_plusPasa/braker.gff3 \
--kill_list red5.from.Russell.cp.09.09.vs.braker3.kill \
--output red5.braker.up.v3.gff3

#now clean up the gff3 - hopefully this will re-adjust the genes to match the transcripts remaining

agat_convert_sp_gxf2gxf.pl -g red5.braker.up.v3.gff3 -o red5.braker.up.clean.v3.gff3

agat_sp_merge_annotations.pl --gff /workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.gff3 \
--gff red5.braker.up.clean.v3.gff3 \
--out braker.liftoff.agat.merge.v3.gff3
EOF

sbatch agat.extract.merge.v3.sl
```

``` r
library(stringr)
library(dplyr)
braker.liftoff.agat.merge.v3.gff3 = read.csv(file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/braker.liftoff.agat.merge.v3.gff3",sep = "\t",skip = 3,header = FALSE)

braker.liftoff.agat.merge.v3.gff3 %>% filter(V3 == "gene") %>% group_by(V2 )  %>% summarise(Count = n())

braker.liftoff.agat.merge.v3.gff3 %>%
  filter(V3 == "gene", V2 == "Liftoff") %>% dplyr::select(V9) %>% 
  mutate(Extracted_ID = str_extract(V9, "ID=[^;]+")) %>%
  mutate(Extracted_ID = str_replace(Extracted_ID, "ID=", ""))  


read.csv(file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.gff3",sep = "\t",skip = 3,header = FALSE) %>% filter(V3 == "gene") %>% group_by(V2 )  %>% summarise(Count = n())
```

``` bash

#get the names liftoff genes that made it through the merge
awk '$2 == "Liftoff" &&
$3 == "gene"'  /workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/braker.liftoff.agat.merge.v3.gff3 | awk '{ split($9, a, /[=;]/); print a[2] }' | sort > liftoff.in.braker.liftoff.agat.merge.v3

#get the names in the original liftoff
awk '$2 == "Liftoff" &&
$3 == "gene"'  /workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.gff3 | awk '{ split($9, a, /[=;]/); print a[2] }' | sort > red5.from.Russell.names
#find the difference between them 
comm -3 red5.from.Russell.names liftoff.in.braker.liftoff.agat.merge.v3 
```

These are the extra genes below responsible for the difference in
number. They were merged by the agat_sp_merge_annotations.pl tool
because all these genes have multiple isoforms. There doesn’t ssem to be
a way to disable this so will use cat to combine the gff3 files.  

    gRUSV2a.032659.1.PC
    gRUSV2a.051718.1.PC
    gRUSV2a.052364.2.PC
    gRUSV2a.052596.2.PC
    RUSV2a.007162.1.PCg427309.431534
    RUSV2a.026153.1.PCg19089532.19091738
    RUSV2a.030344.1.PCg1720849.1726721

Cat the braker gene models and the liftoff together

``` bash
cd $hdir/merge_purge/

cat /workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/red5.from.Russell.cp.09.09.gff3  red5.braker.up.clean.v3.gff3 >  v3.cat.gff3
head red5.braker.up.clean.v3.gff3 
```

Read in the merged file and check the numbers of genes from each source

``` r
v3.cat.gff3 = read.csv(file = "../merge_purge/v3.cat.gff3",sep = "\t",skip = 3,header = FALSE)
v3.cat.gff3 %>%  filter(V3 == "gene") %>% group_by(V2 )  %>% summarise(Count = n())
```

These now match the number of input genes

<table style="width:42%;">
<colgroup>
<col style="width: 25%" />
<col style="width: 16%" />
</colgroup>
<tbody>
<tr class="odd">
<td><strong>Gene source</strong></td>
<td><strong>Count</strong></td>
</tr>
<tr class="even">
<td>AUGUSTUS</td>
<td>5209</td>
</tr>
<tr class="odd">
<td>GeneMark.hmm3</td>
<td>671</td>
</tr>
<tr class="even">
<td>Liftoff</td>
<td>30693</td>
</tr>
<tr class="odd">
<td>gmst</td>
<td>644</td>
</tr>
</tbody>
</table>

I would like like to see if the generic agat_convert_sp_gxf2gxf.pl clean
up tool will do anything annoying like merge isoforms

``` bash

cd $hdir/merge_purge

cat << 'EOF' >  clean.test.sl
#!/bin/bash -e
#SBATCH -J clean.test
#SBATCH --output=clean.test.%j.out
#SBATCH --error=clean.test.%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00


#see if we can re-order the gff nucely without doing any merging

module load conda
conda activate /workspace/hrtjbs/agat
agat_convert_sp_gxf2gxf.pl -g /workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/v3.cat.gff3 -o v3.cat.clean.gff3
EOF
sbatch clean.test.sl
```

``` r
v3.cat.clean.gff3 = read.csv("/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/v3.cat.clean.gff3",header = FALSE,sep = "\t",skip = 3)


v3.cat.clean.gff3 %>%  filter(V3 == "gene") %>% group_by(V2 )  %>% summarise(Count = n())
```

It worked. The breakdown is the same as pre clean - up.

<table style="width:33%;">
<colgroup>
<col style="width: 22%" />
<col style="width: 11%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">Gene source</th>
<th style="text-align: right;">Count</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">AUGUSTUS</td>
<td style="text-align: right;">5209</td>
</tr>
<tr class="even">
<td style="text-align: left;">GeneMark.hmm3</td>
<td style="text-align: right;">671</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Liftoff</td>
<td style="text-align: right;">30693</td>
</tr>
<tr class="even">
<td style="text-align: left;">gmst</td>
<td style="text-align: right;">644</td>
</tr>
</tbody>
</table>

``` r
setwd("/workspace/hrtjbs/pangene_2024/pangene/prototyping/liftoff/")
quarto::quarto_render(input = "Pangene_gene.selection.qmd",output_file = "Pangene_gene.selection.md",execute = FALSE,output_format = "md")
```
