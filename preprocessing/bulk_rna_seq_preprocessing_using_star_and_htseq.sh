#!/bin/bash

# Note that I didn’t use a pipeline for preprocessing—I ran each command manually. So, strictly speaking, this script isn’t a pipeline, but rather a record of the commands I used during preprocessing.

# This shows the preprocessing steps applied to samples SRR21570352 through SRR21570362.

# Set paths to tools and directories (update these accordingly)
SRATOOLKIT_PATH="path/to/sratoolkit.3.1.1-win64/bin"
STAR_PATH="path/to/STAR-2.7.11b/source"
SAMTOOLS_PATH="path/to/samtools"
HTSEQ_PATH="path/to/htseq"
INFER_EXP_PATH="path/to/infer_experiment.py"
GENOME_DIR="path/to/index_for_elife_gencode"
FASTQ_DIR="path/to/fastq_current/RNA_seq_file"
OUTPUT_DIR="path/to/bam_for_elife_gencode"
GTF_FILE="path/to/gencode.v46.primary_assembly.annotation.gtf"
BED_FILE="path/to/gencode.v46.primary_assembly.annotation.bed"

# Define SRR numbers
SRR_LIST=(62 61 60 59 58 57 56 55 54 53 52)

# Download SRR files
echo "Downloading SRR files..."
for i in "${SRR_LIST[@]}"; do
    $SRATOOLKIT_PATH/Prefetch SRR215703$i
done

# Validate download
echo "Validating downloads..."
for i in "${SRR_LIST[@]}"; do
    $SRATOOLKIT_PATH/vdb-validate SRR215703$i
done

# Convert SRR files to FASTQ
echo "Converting SRR files to FASTQ..."
for i in "${SRR_LIST[@]}"; do
    $SRATOOLKIT_PATH/fasterq-dump SRR215703$i
done

# Create STAR index
echo "Creating STAR genome index..."
~/Desktop/STAR/source/STAR --runThreadN 7 --runMode genomeGenerate --genomeDir ~/Desktop/for_STAR_indexing_gencode/index \
    --genomeFastaFiles ~/Desktop/for_STAR_indexing_gencode/GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile ~/Desktop/for_STAR_indexing_gencode/gencode.v46.primary_assembly.annotation.gtf \
    --sjdbOverhang 150

# Mapping using STAR
echo "Mapping reads to reference genome..."
MAPPING_LIST=(21570352 21570353 21570354 21570355 21570356 21570357 21570358 21570359 21570360 21570361 21570362)
for i in "${MAPPING_LIST[@]}"; do
    $STAR_PATH/STAR --outFileNamePrefix $OUTPUT_DIR/SRR${i} --outSAMtype BAM SortedByCoordinate --runThreadN 10 \
        --genomeDir $GENOME_DIR --readFilesIn $FASTQ_DIR/SRR${i}_1.fastq $FASTQ_DIR/SRR${i}_2.fastq
done

# Sorting BAM files by name
echo "Sorting BAM files by name..."
SORT_LIST=(21570352 21570353 21570354 21570355 21570356 21570357 21570358 21570359 21570360 21570361 21570362)
for i in "${SORT_LIST[@]}"; do
    $SAMTOOLS_PATH/samtools sort -n -o $OUTPUT_DIR/SRR${i}_name_sorted.bam $OUTPUT_DIR/SRR${i}.bam
done

# Convert GTF to BED
echo "Converting GTF to BED..."
gff2bed < ~/Desktop/for_STAR_indexing/gencode.v46.primary_assembly.annotation.gtf > $BED_FILE

# Determine strandedness
echo "Determining strandedness using infer_experiment.py..."
for i in "${SORT_LIST[@]}"; do
    python3 $INFER_EXP_PATH -r $BED_FILE -i $OUTPUT_DIR/SRR${i}_name_sorted.bam
done

# You should determine the appropriate strandedness based on the results from the 'Determine strandedness' step above, and apply it using the -s parameter in the command below (the example below uses "no" for the -s parameter).

# Counting reads with HTSeq
echo "Counting reads with HTSeq..."
htseq-count -f bam -r name -s no -t exon -i gene_id --additional-attr=gene_name -m union --nonunique=none \
    -c path/to/output_counts.tsv path/to/input_reads.bam $GTF_FILE

