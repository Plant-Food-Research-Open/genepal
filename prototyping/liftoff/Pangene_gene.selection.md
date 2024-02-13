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
           'donorgff' = '/workspace/hrtjbs/Annotation_Pipeline_2023/liftOff_tests/refseqs/RU01.20221117162301.no.extras.primary.gff3squeue -u hrtjbs
           ',
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
