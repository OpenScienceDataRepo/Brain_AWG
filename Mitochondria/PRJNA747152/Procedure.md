# FASTQ -> Count Matrices

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

## 2. FASTQC -> BAM
I moved the resulting FastQC files to a folder called "FastQC". I then used the [HISAT2 package](https://daehwankimlab.github.io/hisat2/download/) to align the reads to the reference genome.

Note that some packages have to be installed and initialized if you don't have it already:
  [homebrew](https://brew.sh/)
  [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main)
  [hisat2](https://daehwankimlab.github.io/hisat2/download/)
  [samtools](https://www.htslib.org/)

The following bash file was created for simplicity, titled "run_alignment_safe.sh":
```
#!/bin/bash

# Run this script with: bash run_alignment_safe.sh

# Paths
FASTQ_DIR="/Volumes/Marians_SSD/ADBR_Mito/PRJNA747152/FASTQ"
OUTPUT_DIR="/Volumes/Marians_SSD/ADBR_Mito/PRJNA747152"
INDEX="/Volumes/Marians_SSD/dm6/genome"

# Activate conda environment
echo "Activating rnaseq_env..."
# Ensure conda.sh is sourced so 'conda activate' works
source /opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh
conda activate rnaseq_env

# Confirm binaries
echo "Using HISAT2 at: $(which hisat2)"
echo "Using SAMtools at: $(which samtools)"

# Count total FASTQs
TOTAL=$(ls "$FASTQ_DIR"/*.fastq.gz 2>/dev/null | wc -l)
BAM_COUNT=0

# Loop through FASTQ files
for fq in "$FASTQ_DIR"/*.fastq.gz; do
    base=$(basename "$fq" .fastq.gz)
    bam="$OUTPUT_DIR/${base}.bam"

    if [ -f "$bam" ]; then
        echo "[$base] BAM exists, skipping."
        BAM_COUNT=$((BAM_COUNT + 1))
        continue
    fi

    echo "[$base] Aligning FASTQ..."
    hisat2 -x "$INDEX" -U "$fq" -S "$OUTPUT_DIR/${base}.sam"

    echo "[$base] Sorting SAM to BAM..."
    samtools sort -o "$bam" "$OUTPUT_DIR/${base}.sam"
    rm "$OUTPUT_DIR/${base}.sam"

    BAM_COUNT=$((BAM_COUNT + 1))
    echo "=== Progress: $BAM_COUNT / $TOTAL BAMs completed ==="
done

echo "All files processed."
date
```

## Citations
Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4
