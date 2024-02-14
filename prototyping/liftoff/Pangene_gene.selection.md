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

Rename outfile (forget to specify this properly in the script)

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


write.csv(x = merged.summary,file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/merged.gff.summary.csv")
write.csv(x = donor.summary,file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/donor.gff.summary.csv")
write.csv(x = braker.summary,file = "/workspace/hrtjbs/pangene_2024/pangene/prototyping/merge_purge/braker.gff.summary,csv")
```

``` r
quarto::quarto_render(input = "Pangene_gene.selection.qmd",output_file = "Pangene_gene.selection.md",execute = FALSE,output_format = "md")
```
