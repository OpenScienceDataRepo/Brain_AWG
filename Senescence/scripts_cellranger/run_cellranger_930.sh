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
show_progress 0 8 "SRX28491911"
cd /home/easlinger/data/OSD-930
cellranger count --create-bam=false --id=SRX28491911 --sample=SRX28491911 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-930  # SRX28491911
show_progress 1 8 "done: SRX28491911"
show_progress 1 8 "SRX28491912"
cd /home/easlinger/data/OSD-930
cellranger count --create-bam=false --id=SRX28491912 --sample=SRX28491912 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-930  # SRX28491912
show_progress 2 8 "done: SRX28491912"
show_progress 2 8 "SRX28491914"
cd /home/easlinger/data/OSD-930
cellranger count --create-bam=false --id=SRX28491914 --sample=SRX28491914 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-930  # SRX28491914
show_progress 3 8 "done: SRX28491914"
show_progress 3 8 "SRX28491915"
cd /home/easlinger/data/OSD-930
cellranger count --create-bam=false --id=SRX28491915 --sample=SRX28491915 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-930  # SRX28491915
show_progress 4 8 "done: SRX28491915"
show_progress 4 8 "SRX28491917"
cd /home/easlinger/data/OSD-930
cellranger count --create-bam=false --id=SRX28491917 --sample=SRX28491917 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-930  # SRX28491917
show_progress 5 8 "done: SRX28491917"
show_progress 5 8 "SRX28491918"
cd /home/easlinger/data/OSD-930
cellranger count --create-bam=false --id=SRX28491918 --sample=SRX28491918 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-930  # SRX28491918
show_progress 6 8 "done: SRX28491918"
show_progress 6 8 "SRX28491919"
cd /home/easlinger/data/OSD-930
cellranger count --create-bam=false --id=SRX28491919 --sample=SRX28491919 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-930  # SRX28491919
show_progress 7 8 "done: SRX28491919"
show_progress 7 8 "SRX28491920"
cd /home/easlinger/data/OSD-930
cellranger count --create-bam=false --id=SRX28491920 --sample=SRX28491920 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-930  # SRX28491920
show_progress 8 8 "done: SRX28491920"
printf "\n" >&2
