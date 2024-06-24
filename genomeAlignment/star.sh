#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --mem=60GB
#SBATCH --time=12:00:00
#SBATCH --job-name=star-genome-mapping
#SBATCH --output=/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/scripts/genomeAlignment/%j.out
#SBATCH --error=/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/scripts/genomeAlignment/%j.err







module load STAR

STAR --runThreadN 29 --runMode genomeGenerate --genomeDir /home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/starGenome --genomeFastaFiles /home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/genomeCellRanger/refdata-gex-GRCh38-2020-A/refdata-gex-GRCh38-2020-A/fasta/genome.fa --sjdbGTFfile /home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/genomeCellRanger/refdata-gex-GRCh38-2020-A/refdata-gex-GRCh38-2020-A/genes/genes.gtf

