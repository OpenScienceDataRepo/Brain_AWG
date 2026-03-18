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

hisat2 -q --rna-strandness R \
-x /Volumes/Marians\ SSD/dm6/genome \ -U /Volumes/Marians\ SSD/ADBR\ Mito/PRJNA747152/FASTQ/sample.fastq | samtools sort -o /Volumes/Marians\ SSD/ADBR\ Mito/PRJNA747152/sample.bam
```

## Citations
Kim, D., Paggi, J.M., Park, C. et al. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol 37, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4
