#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --partition=orfoz
#SBATCH --output=logs/star_index_%j.out
#SBATCH --error=logs/star_index_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=192G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dinc_s@ibu.edu.tr

set -e

cd /arf/scratch/sdinc/qe_rnaseq_truba

echo "Job started on:"
hostname
date

echo "Current folder:"
pwd

echo "Creating output folder:"
mkdir -p genome_index
ls -lh genome_index

echo "Reference files:"
ls -lh ref

echo "Creating STAR genome index"

apptainer exec /arf/scratch/sdinc/qe_rnaseq_truba/containers/rnaseq_star_tools.sif STAR \
  --runMode genomeGenerate \
  --genomeDir genome_index \
  --genomeFastaFiles ref/complete.genomic.fna \
  --sjdbGTFfile ref/complete.genomic.gtf \
  --sjdbOverhang 149 \
  --genomeSAindexNbases 13 \
  --runThreadN 56

echo "STAR genome index finished"
date

echo "Index files:"
ls -lh genome_index
