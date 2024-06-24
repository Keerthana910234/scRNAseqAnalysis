#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --mem=60GB
#SBATCH --time=12:00:00
#SBATCH --job-name=genome-alignment-to-reference-2020
#SBATCH --output=/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/scripts/genomeAlignment/logs/%x-%j.out
#SBATCH --error=/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/scripts/genomeAlignment/logs/%x-%j.err

module load STAR

# Path to the manifest file
MANIFEST=/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/rawData/bulkRNAseqData/manifest.tsv

# Read each line in the manifest
while IFS=$'\t' read -r col1 col2 col3
do
    # Define output prefix using the sample identifier from col3
    OUTPUT_PREFIX="/home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/scripts/genomeAlignment/outputNew/${col3}_"

    # Ensure output directory exists
    mkdir -p "$(dirname "$OUTPUT_PREFIX")"

    # Run STAR for each sample
    STAR \
    --outSAMattributes All \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --readFilesCommand zcat \
    --runThreadN 30 \
    --sjdbGTFfile /home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/genomeCellRanger/refdata-gex-GRCh38-2020-A/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
    --outReadsUnmapped Fastx \
    --outMultimapperOrder Random \
    --outWigType wiggle \
    --genomeDir /home/mzo5929/Keerthana/subiaHanxiaoDataAnalysis/bulkRNAseqData/starGenome \
    --readFilesIn "$col1" "$col2" \
    --outFileNamePrefix "$OUTPUT_PREFIX"
done < "$MANIFEST"

