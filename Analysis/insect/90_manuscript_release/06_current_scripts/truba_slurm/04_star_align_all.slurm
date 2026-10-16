#!/bin/bash
#SBATCH --job-name=star_align_all
#SBATCH --partition=orfoz
#SBATCH --output=logs/star_align_all_%j.out
#SBATCH --error=logs/star_align_all_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=128G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dinc_s@ibu.edu.tr

set -e

ulimit -n 65535

cd /arf/scratch/sdinc/qe_rnaseq_truba

echo "Job started on:"
hostname
date

echo "Current folder:"
pwd

echo "Open file limit:"
ulimit -n

echo "STAR index folder:"
ls -lh genome_index

echo "Input FASTQ files:"
ls -lh raw_data

echo "Creating alignment output folders"
mkdir -p alignments/L1
mkdir -p alignments/L3
mkdir -p alignments/L4
mkdir -p alignments/P3
mkdir -p alignments/P4
mkdir -p alignments/P5


echo "Starting STAR alignment for L1"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif STAR --genomeDir genome_index --readFilesIn raw_data/L1.R1.fq.gz raw_data/L1.R2.fq.gz --readFilesCommand zcat --runThreadN 28 --outFileNamePrefix alignments/L1/ --outSAMtype BAM Unsorted

echo "Sorting BAM for L1"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools sort -@ 28 -m 3G -o alignments/L1/L1.sorted.bam alignments/L1/Aligned.out.bam

echo "Indexing BAM for L1"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools index alignments/L1/L1.sorted.bam

echo "Removing unsorted BAM for L1"
rm -f alignments/L1/Aligned.out.bam

echo "L1 finished"
date


echo "Starting STAR alignment for L3"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif STAR --genomeDir genome_index --readFilesIn raw_data/L3.R1.fq.gz raw_data/L3.R2.fq.gz --readFilesCommand zcat --runThreadN 28 --outFileNamePrefix alignments/L3/ --outSAMtype BAM Unsorted

echo "Sorting BAM for L3"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools sort -@ 28 -m 3G -o alignments/L3/L3.sorted.bam alignments/L3/Aligned.out.bam

echo "Indexing BAM for L3"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools index alignments/L3/L3.sorted.bam

echo "Removing unsorted BAM for L3"
rm -f alignments/L3/Aligned.out.bam

echo "L3 finished"
date


echo "Starting STAR alignment for L4"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif STAR --genomeDir genome_index --readFilesIn raw_data/L4.R1.fq.gz raw_data/L4.R2.fq.gz --readFilesCommand zcat --runThreadN 28 --outFileNamePrefix alignments/L4/ --outSAMtype BAM Unsorted

echo "Sorting BAM for L4"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools sort -@ 28 -m 3G -o alignments/L4/L4.sorted.bam alignments/L4/Aligned.out.bam

echo "Indexing BAM for L4"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools index alignments/L4/L4.sorted.bam

echo "Removing unsorted BAM for L4"
rm -f alignments/L4/Aligned.out.bam

echo "L4 finished"
date


echo "Starting STAR alignment for P3"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif STAR --genomeDir genome_index --readFilesIn raw_data/P3.R1.fq.gz raw_data/P3.R2.fq.gz --readFilesCommand zcat --runThreadN 28 --outFileNamePrefix alignments/P3/ --outSAMtype BAM Unsorted

echo "Sorting BAM for P3"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools sort -@ 28 -m 3G -o alignments/P3/P3.sorted.bam alignments/P3/Aligned.out.bam

echo "Indexing BAM for P3"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools index alignments/P3/P3.sorted.bam

echo "Removing unsorted BAM for P3"
rm -f alignments/P3/Aligned.out.bam

echo "P3 finished"
date


echo "Starting STAR alignment for P4"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif STAR --genomeDir genome_index --readFilesIn raw_data/P4.R1.fq.gz raw_data/P4.R2.fq.gz --readFilesCommand zcat --runThreadN 28 --outFileNamePrefix alignments/P4/ --outSAMtype BAM Unsorted

echo "Sorting BAM for P4"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools sort -@ 28 -m 3G -o alignments/P4/P4.sorted.bam alignments/P4/Aligned.out.bam

echo "Indexing BAM for P4"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools index alignments/P4/P4.sorted.bam

echo "Removing unsorted BAM for P4"
rm -f alignments/P4/Aligned.out.bam

echo "P4 finished"
date


echo "Starting STAR alignment for P5"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif STAR --genomeDir genome_index --readFilesIn raw_data/P5.R1.fq.gz raw_data/P5.R2.fq.gz --readFilesCommand zcat --runThreadN 28 --outFileNamePrefix alignments/P5/ --outSAMtype BAM Unsorted

echo "Sorting BAM for P5"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools sort -@ 28 -m 3G -o alignments/P5/P5.sorted.bam alignments/P5/Aligned.out.bam

echo "Indexing BAM for P5"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif samtools index alignments/P5/P5.sorted.bam

echo "Removing unsorted BAM for P5"
rm -f alignments/P5/Aligned.out.bam

echo "P5 finished"
date


echo "All STAR alignments finished"

echo "Alignment folders:"
ls -lh alignments

echo "Final sorted BAM files:"
ls -lh alignments/*/*.sorted.bam

echo "BAM index files:"
ls -lh alignments/*/*.sorted.bam.bai

echo "STAR final log files:"
ls -lh alignments/*/Log.final.out

echo "Job finished on:"
date
