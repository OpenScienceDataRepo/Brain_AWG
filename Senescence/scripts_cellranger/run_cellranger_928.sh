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
show_progress 0 8 "SRX28491805"
cd /home/easlinger/data/OSD-928
cellranger count --create-bam=false --id=SRX28491805 --sample=SRX28491805 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-928  # SRX28491805
show_progress 1 8 "done: SRX28491805"
show_progress 1 8 "SRX28491806"
cd /home/easlinger/data/OSD-928
cellranger count --create-bam=false --id=SRX28491806 --sample=SRX28491806 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-928  # SRX28491806
show_progress 2 8 "done: SRX28491806"
show_progress 2 8 "SRX28491807"
cd /home/easlinger/data/OSD-928
cellranger count --create-bam=false --id=SRX28491807 --sample=SRX28491807 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-928  # SRX28491807
show_progress 3 8 "done: SRX28491807"
show_progress 3 8 "SRX28491808"
cd /home/easlinger/data/OSD-928
cellranger count --create-bam=false --id=SRX28491808 --sample=SRX28491808 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-928  # SRX28491808
show_progress 4 8 "done: SRX28491808"
show_progress 4 8 "SRX28491809"
cd /home/easlinger/data/OSD-928
cellranger count --create-bam=false --id=SRX28491809 --sample=SRX28491809 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-928  # SRX28491809
show_progress 5 8 "done: SRX28491809"
show_progress 5 8 "SRX28491810"
cd /home/easlinger/data/OSD-928
cellranger count --create-bam=false --id=SRX28491810 --sample=SRX28491810 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-928  # SRX28491810
show_progress 6 8 "done: SRX28491810"
show_progress 6 8 "SRX28491811"
cd /home/easlinger/data/OSD-928
cellranger count --create-bam=false --id=SRX28491811 --sample=SRX28491811 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-928  # SRX28491811
show_progress 7 8 "done: SRX28491811"
show_progress 7 8 "SRX28491812"
cd /home/easlinger/data/OSD-928
cellranger count --create-bam=false --id=SRX28491812 --sample=SRX28491812 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-928  # SRX28491812
show_progress 8 8 "done: SRX28491812"
printf "\n" >&2
