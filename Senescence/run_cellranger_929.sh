#!/bin/bash
commands=(
    "cellranger count --create-bam=false --id=SRX28491814 --sample=SRX28491814 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491814"
    "cellranger count --create-bam=false --id=SRX28491814 --sample=SRX28491814 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491814"
    "cellranger count --create-bam=false --id=SRX28491815 --sample=SRX28491815 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491815"
    "cellranger count --create-bam=false --id=SRX28491815 --sample=SRX28491815 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491815"
    "cellranger count --create-bam=false --id=SRX28491816 --sample=SRX28491816 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491816"
    "cellranger count --create-bam=false --id=SRX28491816 --sample=SRX28491816 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491816"
    "cellranger count --create-bam=false --id=SRX28491849 --sample=SRX28491849 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491849"
    "cellranger count --create-bam=false --id=SRX28491849 --sample=SRX28491849 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491849"
    "cellranger count --create-bam=false --id=SRX28491850 --sample=SRX28491850 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491850"
    "cellranger count --create-bam=false --id=SRX28491850 --sample=SRX28491850 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491850"
    "cellranger count --create-bam=false --id=SRX28491851 --sample=SRX28491851 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491851"
    "cellranger count --create-bam=false --id=SRX28491851 --sample=SRX28491851 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491851"
    "cellranger count --create-bam=false --id=SRX28491852 --sample=SRX28491852 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491852"
    "cellranger count --create-bam=false --id=SRX28491852 --sample=SRX28491852 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491852"
    "cellranger count --create-bam=false --id=SRX28491853 --sample=SRX28491853 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491853"
    "cellranger count --create-bam=false --id=SRX28491853 --sample=SRX28491853 --localcores=32 --localmem=64 --transcriptome=/home/easlinger/data/reference_genomes/refdata-GRCm39-M36-mouse/mkref_GRCm39 --fastqs=/home/easlinger/data/OSD-929_raw --output-dir=/home/easlinger/data/OSD-929  # SRX28491853"
)

total=${#commands[@]}

draw_progress() {
    local current=$1
    local cols=$(tput cols)
    local bar_width=$((cols - 20))
    local filled=$((current * bar_width / total))
    local bar=$(printf '█%.0s' $(seq 1 $filled))$(printf '░%.0s' $(seq 1 $((bar_width - filled))))
    tput sc
    tput cup $(tput lines) 0
    printf "[%s] %d/%d" "$bar" "$current" "$total"
    tput rc
}

for i in "${!commands[@]}"; do
    draw_progress "$i"
    echo ">>> Running: ${commands[$i]##*#}"
    eval "${commands[$i]%%#*}"
    draw_progress "$((i + 1))"
done

tput cup $(tput lines) 0
printf "\nDone! %d/%d commands completed.\n" "$total" "$total"
