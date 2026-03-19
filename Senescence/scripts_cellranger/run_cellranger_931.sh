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
show_progress 0 8 "SRX28491921"
cd /home/easlinger/data/OSD-931
cellranger count --create-bam=false --id=SRX28491921 --sample=SRX28491921 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-931  # SRX28491921
show_progress 1 8 "done: SRX28491921"
show_progress 1 8 "SRX28491922"
cd /home/easlinger/data/OSD-931
cellranger count --create-bam=false --id=SRX28491922 --sample=SRX28491922 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-931  # SRX28491922
show_progress 2 8 "done: SRX28491922"
show_progress 2 8 "SRX28491971"
cd /home/easlinger/data/OSD-931
cellranger count --create-bam=false --id=SRX28491971 --sample=SRX28491971 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-931  # SRX28491971
show_progress 3 8 "done: SRX28491971"
show_progress 3 8 "SRX28491972"
cd /home/easlinger/data/OSD-931
cellranger count --create-bam=false --id=SRX28491972 --sample=SRX28491972 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-931  # SRX28491972
show_progress 4 8 "done: SRX28491972"
show_progress 4 8 "SRX28491973"
cd /home/easlinger/data/OSD-931
cellranger count --create-bam=false --id=SRX28491973 --sample=SRX28491973 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-931  # SRX28491973
show_progress 5 8 "done: SRX28491973"
show_progress 5 8 "SRX28491974"
cd /home/easlinger/data/OSD-931
cellranger count --create-bam=false --id=SRX28491974 --sample=SRX28491974 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-931  # SRX28491974
show_progress 6 8 "done: SRX28491974"
show_progress 6 8 "SRX28491976"
cd /home/easlinger/data/OSD-931
cellranger count --create-bam=false --id=SRX28491976 --sample=SRX28491976 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-931  # SRX28491976
show_progress 7 8 "done: SRX28491976"
show_progress 7 8 "SRX28491977"
cd /home/easlinger/data/OSD-931
cellranger count --create-bam=false --id=SRX28491977 --sample=SRX28491977 --localcores=16 --localmem=45 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-931  # SRX28491977
show_progress 8 8 "done: SRX28491977"
printf "\n" >&2
