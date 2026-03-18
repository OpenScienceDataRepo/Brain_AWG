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
show_progress 0 8 "SRX28491814"
cd /home/easlinger/data/OSD-929
cellranger count --create-bam=false --id=SRX28491814 --sample=SRX28491814 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929  # SRX28491814
show_progress 1 8 "done: SRX28491814"
show_progress 1 8 "SRX28491815"
cd /home/easlinger/data/OSD-929
cellranger count --create-bam=false --id=SRX28491815 --sample=SRX28491815 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929  # SRX28491815
show_progress 2 8 "done: SRX28491815"
show_progress 2 8 "SRX28491816"
cd /home/easlinger/data/OSD-929
cellranger count --create-bam=false --id=SRX28491816 --sample=SRX28491816 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929  # SRX28491816
show_progress 3 8 "done: SRX28491816"
show_progress 3 8 "SRX28491849"
cd /home/easlinger/data/OSD-929
cellranger count --create-bam=false --id=SRX28491849 --sample=SRX28491849 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929  # SRX28491849
show_progress 4 8 "done: SRX28491849"
show_progress 4 8 "SRX28491850"
cd /home/easlinger/data/OSD-929
cellranger count --create-bam=false --id=SRX28491850 --sample=SRX28491850 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929  # SRX28491850
show_progress 5 8 "done: SRX28491850"
show_progress 5 8 "SRX28491851"
cd /home/easlinger/data/OSD-929
cellranger count --create-bam=false --id=SRX28491851 --sample=SRX28491851 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929  # SRX28491851
show_progress 6 8 "done: SRX28491851"
show_progress 6 8 "SRX28491852"
cd /home/easlinger/data/OSD-929
cellranger count --create-bam=false --id=SRX28491852 --sample=SRX28491852 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929  # SRX28491852
show_progress 7 8 "done: SRX28491852"
show_progress 7 8 "SRX28491853"
cd /home/easlinger/data/OSD-929
cellranger count --create-bam=false --id=SRX28491853 --sample=SRX28491853 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929  # SRX28491853
show_progress 8 8 "done: SRX28491853"
printf "\n" >&2
