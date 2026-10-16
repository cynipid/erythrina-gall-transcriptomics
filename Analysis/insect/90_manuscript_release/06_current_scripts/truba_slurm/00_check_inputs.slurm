#!/bin/bash
#SBATCH --job-name=check_inputs
#SBATCH --partition=debug
#SBATCH --output=logs/build_rnaseq_sif_check_inputs_%j.out
#SBATCH --error=logs/build_rnaseq_sif_check_inputs.%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=56
#SBATCH --mem=16G

cd /arf/scratch/sdinc/qe_rnaseq_truba/

echo "Job started on:"
hostname
date

echo "Current folder:"
pwd

echo "Raw data files:"
ls -lh raw_data

echo "Reference files:"
ls -lh ref

echo "Testing gzip files"

apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/L1.R1.fq.gz
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/L1.R2.fq.gz

apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/L3.R1.fq.gz
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/L3.R2.fq.gz

apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/L4.R1.fq.gz
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/L4.R2.fq.gz

apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/P3.R1.fq.gz
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/P3.R2.fq.gz

apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/P4.R1.fq.gz
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/P4.R2.fq.gz

apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/P5.R1.fq.gz
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif pigz -t -p 56 raw_data/P5.R2.fq.gz

echo "All gzip tests finished"
date
