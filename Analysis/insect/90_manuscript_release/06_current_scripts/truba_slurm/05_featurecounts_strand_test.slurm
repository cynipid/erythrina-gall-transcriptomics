#!/bin/bash
#SBATCH --job-name=fc_strand_test
#SBATCH --partition=orfoz
#SBATCH --output=logs/featurecounts_strand_test_%j.out
#SBATCH --error=logs/featurecounts_strand_test_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=64G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dinc_s@ibu.edu.tr

set -e

cd /arf/scratch/sdinc/qe_rnaseq_truba

echo "Job started on:"
hostname
date

echo "Current folder:"
pwd

echo "Creating counts folder"
mkdir -p counts

echo "Input BAM files:"
ls -lh alignments/*/*.sorted.bam

echo "Annotation file:"
ls -lh ref/complete.genomic.gff

echo "Running featureCounts unstranded test: -s 0"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif featureCounts -T 28 -p --countReadPairs -B -C -F GFF -t exon -g locus_tag -s 0 -a ref/complete.genomic.gff -o counts/gene_counts_s0.txt alignments/L1/L1.sorted.bam alignments/L3/L3.sorted.bam alignments/L4/L4.sorted.bam alignments/P3/P3.sorted.bam alignments/P4/P4.sorted.bam alignments/P5/P5.sorted.bam

echo "Running featureCounts forward-stranded test: -s 1"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif featureCounts -T 28 -p --countReadPairs -B -C -F GFF -t exon -g locus_tag -s 1 -a ref/complete.genomic.gff -o counts/gene_counts_s1.txt alignments/L1/L1.sorted.bam alignments/L3/L3.sorted.bam alignments/L4/L4.sorted.bam alignments/P3/P3.sorted.bam alignments/P4/P4.sorted.bam alignments/P5/P5.sorted.bam

echo "Running featureCounts reverse-stranded test: -s 2"
apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif featureCounts -T 28 -p --countReadPairs -B -C -F GFF -t exon -g locus_tag -s 2 -a ref/complete.genomic.gff -o counts/gene_counts_s2.txt alignments/L1/L1.sorted.bam alignments/L3/L3.sorted.bam alignments/L4/L4.sorted.bam alignments/P3/P3.sorted.bam alignments/P4/P4.sorted.bam alignments/P5/P5.sorted.bam

echo "featureCounts strand test finished"
date

echo "Summary files:"
ls -lh counts/gene_counts_s*.summary

echo "Assigned counts summary:"
grep "Assigned" counts/gene_counts_s0.txt.summary
grep "Assigned" counts/gene_counts_s1.txt.summary
grep "Assigned" counts/gene_counts_s2.txt.summary

echo "Job finished on:"
date
