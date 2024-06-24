#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --mem=60GB
#SBATCH --time=12:00:00
#SBATCH --job-name=genome-alignment-to-reference
#SBATCH --output=/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/scripts/genomeAlignment/logs/%X-%j.out
#SBATCH --error=/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/scripts/genomeAlignment/logs/%X-%j.err





module load STAR

STAR \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--readFilesCommand zcat \
--runThreadN 30 \
--sjdbGTFfile /home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/genome/gencode.v46.primary_assembly.annotation.gtf \
--outReadsUnmapped Fastx \
--outMultimapperOrder Random \
--outWigType wiggle \
--genomeDir /home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/starGenome \
--readFilesManifest /home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/rawData/bulkRNAseqData/manifest.tsv                                                                                 
