#!/bin/bash
#SBATCH --job-name=build_rnaseq_sif
#SBATCH --partition=debug
#SBATCH --output=logs/build_rnaseq_sif_%j.out
#SBATCH --error=logs/build_rnaseq_sif_%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --mem=64G

cd /arf/scratch/sdinc/qe_rnaseq_truba/containers

echo "Job started on:"
hostname
date

echo "Building final SIF container..."
apptainer build --force rnaseq_star_tools.sif rnaseq_miniconda_sandbox

echo "Build finished:"
date

ls -lh rnaseq_star_tools.sif
