#!/bin/bash -e


#SBATCH --job-name PANGENE
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output pangene_pfr.stdout
#SBATCH --error pangene_pfr.stderr
#SBATCH --mem=4G

ml apptainer/1.1
ml nextflow/23.04.4

export APPTAINER_BINDPATH="$APPTAINER_BINDPATH,/workspace/$USER/tmp:/tmp"
export TMPDIR="/workspace/$USER/tmp"

nextflow main.nf -profile slurm -resume