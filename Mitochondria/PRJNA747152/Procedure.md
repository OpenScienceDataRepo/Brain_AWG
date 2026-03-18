# FASTQ -> FASTQC -> HISAT2 -> Count Matrices

## 1. FASTQ -> FASTQC

The FASTQ files for this dataset were downloaded from [this directory](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA747152). They were then converted to Count Matrices using the [FastQC package](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

To install the FastQC Package, I used the following code:
```
# changing directory to where my installation is
cd /Volumes/Marians\ SSD

# unzipping and moving the installation
unzip fastqc_v0.12.1.zip
mv FastQC ~/FastQC
cd ~/FastQC

# make fastqc executable
chmod 755 fastqc

# testing fastqc (if it launches, then it's ready to go)
./fastqc

# make FastQC available globally (using symlink)
sudo ln -s /Users/mariannavarro/FastQC/fastqc /usr/local/bin/fastqc

# test
fastqc --help
```

Once that was set up, I ran fastqc on the PRJNA747152 FastQ files using the following code:
```
fastqc /Volumes/Marians\ SSD/ADBR\ Mito/PRJNA747152/FASTQ/* -o /Volumes/Marians\ SSD/ADBR\ Mito/PRJNA747152
```

## 2. FASTQC -> HISAT2
I moved the resulting FastQC files to a folder called "FastQC". I then used the [HISAT2 package](https://daehwankimlab.github.io/hisat2/download/) to align the reads to the reference genome.

Note that some packages have to be installed and initialized if you don't have it already:
  [homebrew](https://brew.sh/)
  [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)
  [hisat2](https://daehwankimlab.github.io/hisat2/download/)
  [samtools](https://www.htslib.org/)

The following bash file was created for simplicity, titled "run_alignment_with_monitor.sh":
```
#!/bin/bash

# RNA-seq alignment monitor

# Check that we are running under the correct environment
echo "Running in Conda environment: $CONDA_DEFAULT_ENV"
echo "Using HISAT2 at: $(which hisat2)"
echo "Using SAMtools at: $(which samtools)"
echo ""

# Directories (quote paths to handle spaces)
FASTQ_DIR="$HOME/marians_ssd/ADBR_Mito/PRJNA747152/FASTQ"
OUTPUT_DIR="$HOME/marians_ssd/ADBR_Mito/PRJNA747152"

# Count total FASTQ files
TOTAL=$(ls "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | wc -l)
BAM_COUNT=$(ls "$OUTPUT_DIR"/*.bam 2>/dev/null | wc -l)
SAM_COUNT=$(ls "$OUTPUT_DIR"/*.sam 2>/dev/null | wc -l)

# List FASTQs that still need alignment
echo "Files still needing alignment/BAM:"
for fq in "$FASTQ_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    if [ ! -f "$OUTPUT_DIR/${base}.bam" ]; then
        echo "  $base"
    fi
done
echo ""

# Progress summary
echo "=== Alignment Progress ==="
echo "Total FASTQ files: $TOTAL"
echo "BAM files completed: $BAM_COUNT"
echo "SAM files present (not yet BAM): $SAM_COUNT"
echo "Remaining files: $((TOTAL - BAM_COUNT))"
echo "=========================="
echo ""

# Start alignment
echo "Starting alignment for $TOTAL FASTQ files..."
echo ""

for fq in "$FASTQ_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    bam_file="$OUTPUT_DIR/${base}.bam"

    if [ -f "$bam_file" ]; then
        echo "[$base] BAM exists, skipping."
        continue
    fi

    echo "[$base] Aligning FASTQ..."
    hisat2 -x "$OUTPUT_DIR/genome_index" -U "$fq" -S "$OUTPUT_DIR/${base}.sam"
    
    if [ $? -ne 0 ]; then
        echo "Error: HISAT2 failed for $base. Skipping..."
        continue
    fi

    echo "[$base] Sorting SAM to BAM..."
    samtools sort -o "$bam_file" "$OUTPUT_DIR/${base}.sam"
    
    if [ $? -ne 0 ]; then
        echo "Error: SAMtools sort failed for $base."
        continue
    fi

    echo "[$base] Done."
    echo ""
done

# Final summary
BAM_COUNT=$(ls "$OUTPUT_DIR"/*.bam 2>/dev/null | wc -l)
echo "All files processed."
echo "Total BAM files: $BAM_COUNT"
date
```
Run it in the terminal with: ```conda run -n rnaseq_env "/Volumes/Marians SSD/run_alignment_with_monitor.sh"```
## Citations
Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4
