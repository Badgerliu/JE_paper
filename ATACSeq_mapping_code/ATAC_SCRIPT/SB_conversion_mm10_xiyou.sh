#!/usr/bin/env bash
# Add parameter handling
while getopts ":g:o:" opt; do
  case $opt in
    g) GENOME_FA="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    *) echo "Usage: $0 -g genome.fa -o output_dir" >&2
       exit 1 ;;
  esac
done

for sample in *.sam; do
    describer=$(basename "${sample}" .sam)
    echo "Processing $describer"
    
    # Convert SAM to BAM
    samtools view -bT "$GENOME_FA" "$sample" > "${describer}.uns.bam"
    
    # Sort and index
    samtools sort -o "${OUTDIR}/${describer}.bam" "${describer}.uns.bam"
    samtools index "${OUTDIR}/${describer}.bam"
    
    # Cleanup
    rm "${describer}.uns.bam"
done 
