#!/bin/bash
# Author Huan Liu
# Author Organization Wuhan University School of Stomatology
# Usage: Working with ATAC-seq raw fq.gz files and resulted in bam, bw and peak files for whole peaks and NFRs, EXCLUDING MNR, DNR etc.
# Date: 2022-10-02

set -euo pipefail

# Configuration Section - Users should modify these paths
# -------------------------------------------------------
REF_GENOME_DIR="${REF_GENOME_DIR:-/path/to/your/refergenome/mm10_scATAC/fasta}"  # Set default or use environment variable
TRIM_JAR="${TRIM_JAR:-/path/to/Trimmomatic/trimmomatic-0.39.jar}"  # Trimmomatic JAR
ADAPTERS="${ADAPTERS:-/path/to/Trimmomatic/adapters/NexteraPE-PE.fa}"  # Adapter sequences
ATAC_SCRIPT_DIR="${ATAC_SCRIPT_DIR:-/path/to/ATAC_scripts}"  # Contains helper scripts
CONDA_ENV="${CONDA_ENV:-ATAC}"  # Name of conda environment
THREADS="${THREADS:-16}"  # Default number of threads
BOWTIE2_INDEX="${BOWTIE2_INDEX:-$REF_GENOME_DIR/genome}"  # Bowtie2 index path
GENOME_FA="${GENOME_FA:-$REF_GENOME_DIR/genome.fa}"  # Reference genome FASTA
MARKDUPS_JAR="${MARKDUPS_JAR:-$CONDA_PREFIX/share/picard-2.27.5-0/lib/picard.jar}"  # Auto-detect from conda
MACROPEAK_SIZE="${MACROPEAK_SIZE:-2652783500}"  # mm10 effective genome size
NFRPEAK_SIZE="${NFRPEAK_SIZE:-1.87e9}"          # NFR effective genome size

# Directory Structure Setup
OUTPUT_DIR="./analysis"
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# Initialize logging
LOG_FILE="$LOG_DIR/ATAC_pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Runtime Information
echo "### ATAC-seq Pipeline v1.2 ###"
echo "Started: $(date)"
echo "Conda Environment: $CONDA_ENV"
echo "Reference Genome: $GENOME_FA"
echo "Bowtie2 Index: $BOWTIE2_INDEX"
echo "Output Directory: $OUTPUT_DIR"

# Validate required commands
check_command() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "Error: $1 not found!"
        exit 1
    }
}
check_command java
check_command bowtie2
check_command samtools
check_command picard

# Activate Conda Environment
eval "$(conda shell.bash hook)"
if conda activate "$CONDA_ENV"; then
    echo "Activated conda environment: $CONDA_ENV"
else
    echo "Failed to activate conda environment: $CONDA_ENV"
    exit 1
fi

# Main Pipeline
# Part 1: Trimming with Trimmomatic
echo "[1/5] Trimming with Trimmomatic"
mkdir -p "$OUTPUT_DIR/trimmed/paired" "$OUTPUT_DIR/trimmed/unpaired"

for i in *.fq.gz; do
    base=$(basename "$i" | rev | cut -c 10- | rev)
    echo "Processing $base..."
    
    java -jar "$TRIM_JAR" PE -phred33 \
        "${base}_R1.fq.gz" "${base}_R2.fq.gz" \
        "$OUTPUT_DIR/trimmed/paired/${base}_R1_paired.fq.gz" \
        "$OUTPUT_DIR/trimmed/unpaired/${base}_R1_unpaired.fq.gz" \
        "$OUTPUT_DIR/trimmed/paired/${base}_R2_paired.fq.gz" \
        "$OUTPUT_DIR/trimmed/unpaired/${base}_R2_unpaired.fq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10:8:TRUE \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:7
done

# Part 2: Mapping with Bowtie2 (Updated full processing)
echo "[2/5] Mapping with Bowtie2"
mkdir -p "$OUTPUT_DIR/mapped"

for i in "$OUTPUT_DIR/trimmed/paired/"*_paired.fq.gz; do
    base=$(basename "$i" | rev | cut -c 17- | rev)
    echo "Mapping $base..."
    
    bowtie2 -x "$BOWTIE2_INDEX" \
        -1 "$OUTPUT_DIR/trimmed/paired/${base}_R1_paired.fq.gz" \
        -2 "$OUTPUT_DIR/trimmed/paired/${base}_R2_paired.fq.gz" \
        -X 1500 -p "$THREADS" -S "$OUTPUT_DIR/mapped/${base}.sam"
    
    # Filter chromosomes and convert to BAM
    sed '/chrM/d;/random/d;/chrUn/d' "$OUTPUT_DIR/mapped/${base}.sam" > "$OUTPUT_DIR/mapped/${base}_m.sam"
    samtools view -bT "$GENOME_FA" "$OUTPUT_DIR/mapped/${base}_m.sam" > "$OUTPUT_DIR/mapped/${base}_m.bam"
    rm "$OUTPUT_DIR/mapped/${base}.sam" "$OUTPUT_DIR/mapped/${base}_m.sam"
    
    # Sort and mark duplicates
    samtools sort -o "$OUTPUT_DIR/mapped/${base}.bam" "$OUTPUT_DIR/mapped/${base}_m.bam"
    java -jar "$MARKDUPS_JAR" MarkDuplicates \
        I="$OUTPUT_DIR/mapped/${base}.bam" \
        O="$OUTPUT_DIR/mapped/${base}_m_d.bam" \
        M="$OUTPUT_DIR/mapped/${base}_dups.txt" \
        REMOVE_DUPLICATES=true
    samtools index "$OUTPUT_DIR/mapped/${base}_m_d.bam"
    rm "$OUTPUT_DIR/mapped/${base}.bam" "$OUTPUT_DIR/mapped/${base}_m.bam"
done

# Part 3: QC and Insert Size Metrics
echo "[3/5] Quality Control"
mkdir -p "$OUTPUT_DIR/insert_size_distribution"

# Run insert size metrics
"$ATAC_SCRIPT_DIR/insert_xiyou.sh" -d "$OUTPUT_DIR/mapped" -o "$OUTPUT_DIR/insert_size_distribution"

# Part 4: BigWig Generation
echo "[4/5] BigWig Generation"
mkdir -p "$OUTPUT_DIR/bw"

# Run bigwig generation
"$ATAC_SCRIPT_DIR/generate_bw_mm10_xiyou.sh" \
    -d "$OUTPUT_DIR/mapped" \
    -o "$OUTPUT_DIR/bw" \
    -g "$MACROPEAK_SIZE"

# Part 5: Peak Calling
echo "[5/5] Peak Calling"
mkdir -p "$OUTPUT_DIR/peaks"

# Run peak calling
"$ATAC_SCRIPT_DIR/shift_and_call_peak_mm10_xiyou.sh" \
    -d "$OUTPUT_DIR/mapped" \
    -o "$OUTPUT_DIR/peaks" \
    -g "$NFRPEAK_SIZE"

# Final Cleanup and Report
echo "[5/5] Pipeline completed successfully!"
echo "Total runtime: $SECONDS seconds"
echo "Output directory: $OUTPUT_DIR"
echo "Log file: $LOG_FILE" 
