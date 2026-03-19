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
show_progress 0 7 "SRX28491902"
cd /home/easlinger/data/OSD-927
cellranger count --create-bam=false --id=SRX28491902 --sample=SRX28491902 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-927  # SRX28491902
show_progress 1 7 "done: SRX28491902"
show_progress 1 7 "SRX28491903"
cd /home/easlinger/data/OSD-927
cellranger count --create-bam=false --id=SRX28491903 --sample=SRX28491903 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-927  # SRX28491903
show_progress 2 7 "done: SRX28491903"
show_progress 2 7 "SRX28491904"
cd /home/easlinger/data/OSD-927
cellranger count --create-bam=false --id=SRX28491904 --sample=SRX28491904 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-927  # SRX28491904
show_progress 3 7 "done: SRX28491904"
show_progress 3 7 "SRX28491905"
cd /home/easlinger/data/OSD-927
cellranger count --create-bam=false --id=SRX28491905 --sample=SRX28491905 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-927  # SRX28491905
show_progress 4 7 "done: SRX28491905"
show_progress 4 7 "SRX28491906"
cd /home/easlinger/data/OSD-927
cellranger count --create-bam=false --id=SRX28491906 --sample=SRX28491906 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-927  # SRX28491906
show_progress 5 7 "done: SRX28491906"
show_progress 5 7 "SRX28491908"
cd /home/easlinger/data/OSD-927
cellranger count --create-bam=false --id=SRX28491908 --sample=SRX28491908 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-927  # SRX28491908
show_progress 6 7 "done: SRX28491908"
show_progress 6 7 "SRX28491909"
cd /home/easlinger/data/OSD-927
cellranger count --create-bam=false --id=SRX28491909 --sample=SRX28491909 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-927  # SRX28491909
show_progress 7 7 "done: SRX28491909"
printf "\n" >&2
