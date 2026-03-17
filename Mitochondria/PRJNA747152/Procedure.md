# Procedure
## FASTQ -> Count Matrices
The FASTQ files for this dataset were downloaded from [this directory](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA747152). They were then converted to Count Matrices using the [fastqc package](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

```
cd /Volumes/Marians\ SSD     # changing the directory to where my pipeline is
ls                           # verifying this is the correct folder
chmod 755 PRJNA747152.sh     # making the pipeline executible in bash
./PRJNA747152.sh             # running the pipeline
```
