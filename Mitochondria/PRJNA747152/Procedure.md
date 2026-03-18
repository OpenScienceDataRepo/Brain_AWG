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

```
# install homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# install miniconda and initialize
brew install --cask miniconda
conda init zsh
source ~/.zshrc

# create a clean RNA-seq environment and activate it
conda create -n rnaseq_env -c bioconda -c conda-forge hisat2 samtools
conda activate rnaseq_env

# verify installation
hisat2 --version
samtools --version

# create clean aliases
ln -s "/Volumes/Marians SSD" ~/marians_ssd

# run alignment
for file in ~/marians_ssd/ADBR_Mito/PRJNA747152/FASTQ/*.fastq.gz
do
  base=$(basename "$file" .fastq.gz)
  bamfile=~/marians_ssd/ADBR_Mito/PRJNA747152/${base}.bam
  samfile=~/marians_ssd/ADBR_Mito/PRJNA747152/${base}.sam

  # If BAM already exists, skip
  if [ -f "$bamfile" ]; then
      echo "$bamfile already exists, skipping."
      continue
  fi

  # If SAM exists but BAM does not, just convert SAM to BAM
  if [ -f "$samfile" ]; then
      echo "Converting existing $samfile to BAM"
      samtools sort -@ 8 -o "$bamfile" "$samfile"
      continue
  fi

  # Otherwise, run HISAT2 alignment
  echo "Aligning $file with HISAT2"
  hisat2 -p 8 -q --rna-strandness R \
    -x ~/marians_ssd/dm6/genome \
    -U "$file" \
    -S "$samfile"

  # Convert the newly created SAM to BAM
  samtools sort -@ 8 -o "$bamfile" "$samfile"
done

```

## Citations
Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4
