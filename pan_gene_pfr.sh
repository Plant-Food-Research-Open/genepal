#!/bin/bash -e


#SBATCH --job-name PAN_GENE
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output pan_gene_pfr.stdout
#SBATCH --error pan_gene_pfr.stderr
#SBATCH --mem=4G

ml apptainer/1.1
ml nextflow/22.10.4

export TMPDIR="/workspace/$USER/tmp"

nextflow main.nf -profile slurm -resume