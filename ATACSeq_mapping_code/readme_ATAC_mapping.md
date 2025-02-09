# ATAC-seq Processing Pipeline v1.2

A modular pipeline for processing ATAC-seq data from raw FASTQ files to peak calling, optimized for portability and reproducibility.

## Table of Contents
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Pipeline Steps](#pipeline-steps)
- [Output Structure](#output-structure)
- [Customization](#customization)
- [Troubleshooting](#troubleshooting)
- [Helper Scripts](#helper-scripts)

## Requirements <a name="requirements"></a>
- **Conda** (Miniconda or Anaconda)
- **Java Runtime** (v11+)
- **Unix-like system** (Linux/macOS)
- **Minimum Hardware**:
  - 16 CPU cores
  - 32GB RAM
  - 100GB+ storage (depends on dataset size)

## Quick Start <a name="quick-start"></a>
1. Set up environment:
```bash
conda env create -n ATAC -f atac_environment.yaml
conda activate ATAC
```

2. Configure paths:
```bash
export REF_GENOME_DIR="/path/to/mm10/fasta"
export ATAC_SCRIPT_DIR="/path/to/helper/scripts"
```

3. Run pipeline:
```bash
./ATAC_oneClick_mm10_xiyou_v1.2.sh
```

## Configuration <a name="configuration"></a>
### Essential Environment Variables
| Variable | Description | Default |
|----------|-------------|---------|
| `REF_GENOME_DIR` | Path to reference genome directory | *Required* |
| `TRIM_JAR` | Path to Trimmomatic JAR file | `$CONDA_PREFIX/share/trimmomatic.jar` |
| `ATAC_SCRIPT_DIR` | Directory containing helper scripts | *Required* |
| `THREADS` | Number of CPU threads | 16 |

### Genome Configuration
```bash
# mm10 (default)
export MACROPEAK_SIZE="2652783500"
export NFRPEAK_SIZE="1.87e9"

# hg38 example
export MACROPEAK_SIZE="2913022398" 
export NFRPEAK_SIZE="2.7e9"
```

## Pipeline Steps <a name="pipeline-steps"></a>
1. **Trimming**:  
   ```bash
   java -jar $TRIM_JAR PE -phred33 input_R1.fq.gz input_R2.fq.gz ...
   ```

2. **Mapping**:  
   ```bash
   bowtie2 -x $BOWTIE2_INDEX -1 input_R1_paired.fq.gz -2 input_R2_paired.fq.gz ...
   ```

3. **QC Metrics**:
   - Insert size distribution
   - Duplication metrics

4. **Signal Visualization**:
   - BigWig file generation

5. **Peak Calling**:
   - Whole genome peaks
   - NFR (Nucleosome Free Region) peaks

## Output Structure <a name="output-structure"></a>
```
analysis/
├── logs/                   # Detailed processing logs
├── trimmed/                # Trimmed FASTQs
│   ├── paired/             # High-quality paired-end reads
│   └── unpaired/           # Unpaired reads
├── mapped/                 # Alignment results
│   ├── *.bam               # Processed BAM files
│   └── *_dups.txt          # Duplication metrics
├── insert_size_distribution/ # QC plots and metrics
├── bw/                     # Normalized BigWig files
└── peaks/                  # Peak calling results
    ├── narrowPeak/         # MACS2 peak files
    └── summits/            # Peak summit locations
```

## Customization <a name="customization"></a>
### Adapter Sequences
```bash
# For different library prep kits
export ADAPTERS="/path/to/nextera_adapters.fa"
```

### Advanced Options
```bash
# Increase Java memory for large datasets
export JAVA_OPTS="-Xmx32g"

# Custom genome processing
export CHR_FILTER="chrM|random|chrUn"  # Chromosomes to exclude
```

## Troubleshooting <a name="troubleshooting"></a>
### Common Errors
1. **Missing Dependencies**:
   ```bash
   # Verify conda environment
   conda list --name ATAC
   ```

2. **File Not Found Errors**:
   ```bash
   # Validate reference genome paths
   ls $REF_GENOME_DIR/genome.*
   ```

3. **Memory Issues**:
   ```bash
   # Increase memory allocation
   export JAVA_OPTS="-Xmx64g"
   ```

### Quality Control Checks
1. **Trimming Success**:
   - Verify >90% read retention in `trimmed/paired/`

2. **Mapping Rate**:
   - Check Bowtie2 logs for >70% alignment rate

3. **Insert Sizes**:
   - Validate nucleosomal pattern in `insert_size_distribution/`

## Helper Scripts <a name="helper-scripts"></a>
| Script | Usage | Description |
|--------|-------|-------------|
| `insert_xiyou.sh` | `-d input_dir -o output_dir` | Calculate insert size metrics |
| `shift_and_call_peak_mm10_xiyou.sh` | `-d bam_dir -g genome_size` | Perform shift-and-call peak calling |
| `SB_conversion_mm10_xiyou.sh` | `-g genome.fa -o output_dir` | SAM to BAM conversion |

## Support
For assistance, contact:  
[Your Name]  
[Your Email]  
[Issue Tracker Link]
```
