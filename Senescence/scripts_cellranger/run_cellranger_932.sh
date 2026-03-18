#!/bin/bash

show_progress() {
    local current=$1 total=$2 label=$3
    local cols=$(tput cols)
    local bar_width=$((cols - 25))
    local filled=$((current * bar_width / total))
    local empty=$((bar_width - filled))
    local bar=$(printf '█%.0s' $(seq 1 $filled))
    local rest=$(printf '░%.0s' $(seq 1 $empty))
    printf "\r[%s%s] %d/%d %s" "$bar" "$rest" "$current" "$total" "$label" >&2
}
trap 'rm -f __*.mro' EXIT
show_progress 0 8 "SRX28491854"
cd /home/easlinger/data/OSD-932
cellranger count --create-bam=false --id=SRX28491854 --sample=SRX28491854 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-932  # SRX28491854
show_progress 1 8 "done: SRX28491854"
show_progress 1 8 "SRX28491855"
cd /home/easlinger/data/OSD-932
cellranger count --create-bam=false --id=SRX28491855 --sample=SRX28491855 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-932  # SRX28491855
show_progress 2 8 "done: SRX28491855"
show_progress 2 8 "SRX28491857"
cd /home/easlinger/data/OSD-932
cellranger count --create-bam=false --id=SRX28491857 --sample=SRX28491857 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-932  # SRX28491857
show_progress 3 8 "done: SRX28491857"
show_progress 3 8 "SRX28491858"
cd /home/easlinger/data/OSD-932
cellranger count --create-bam=false --id=SRX28491858 --sample=SRX28491858 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-932  # SRX28491858
show_progress 4 8 "done: SRX28491858"
show_progress 4 8 "SRX28491860"
cd /home/easlinger/data/OSD-932
cellranger count --create-bam=false --id=SRX28491860 --sample=SRX28491860 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-932  # SRX28491860
show_progress 5 8 "done: SRX28491860"
show_progress 5 8 "SRX28491861"
cd /home/easlinger/data/OSD-932
cellranger count --create-bam=false --id=SRX28491861 --sample=SRX28491861 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-932  # SRX28491861
show_progress 6 8 "done: SRX28491861"
show_progress 6 8 "SRX28491862"
cd /home/easlinger/data/OSD-932
cellranger count --create-bam=false --id=SRX28491862 --sample=SRX28491862 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-932  # SRX28491862
show_progress 7 8 "done: SRX28491862"
show_progress 7 8 "SRX28491910"
cd /home/easlinger/data/OSD-932
cellranger count --create-bam=false --id=SRX28491910 --sample=SRX28491910 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-932  # SRX28491910
show_progress 8 8 "done: SRX28491910"
printf "\n" >&2
