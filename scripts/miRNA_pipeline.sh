#!/bin/bash

usage() {
    echo "Usage: $0 <fastq_directory> <output_directory> <bwa_index>"
    echo " fastq_directory: Directory containing FASTQ files to analyze"
    echo " output_directory: Output directory for results"
    echo " bwa_index: BWA index prefix (e.g. /path/to/mature_hsa)"
    exit 1
}

if [ $# -lt 3 ]; then
    echo "Error: Missing required arguments"
    usage
fi

fastq_dir="$1"
output_dir="$2"
bwa_index="$3"

if [ ! -d "$fastq_dir" ]; then
    echo "Error: Directory '$fastq_dir' does not exist"
    exit 1
fi

if [ -z "$output_dir" ]; then
    echo "Error: No output directory specified"
    usage
fi

if [ -z "$(find "$fastq_dir" -name "*.fastq" -o -name "*.fq" 2>/dev/null)" ]; then
    echo "Error: No FASTQ files found in directory '$fastq_dir'"
    exit 1
fi

mkdir -p "$output_dir"
mkdir -p "${output_dir}/stats"
mkdir -p "${output_dir}/logs"

process_fastq() {
    local file="$1"
    local output_dir="$2"
    local bwa_index="$3"
    
    file_name=$(basename "${file%.*}")

    intermediate_dir="${output_dir}/intermediate/${file_name}"
    log_dir="${output_dir}/logs"
    mkdir -p "$intermediate_dir"

    umi_extracted_file="${intermediate_dir}/${file_name}_umi_extracted.fastq"
    umi_log="${intermediate_dir}/${file_name}_umi_extraction.log"
    trimmed_file="${intermediate_dir}/${file_name}_trimmed.fastq"
    bwa_output="${intermediate_dir}/${file_name}_bwa.sai"
    bwa_sam="${intermediate_dir}/${file_name}_bwa.sam"
    bwa_bam="${intermediate_dir}/${file_name}_bwa.bam"
    dedup_bam="${intermediate_dir}/${file_name}_deduplicated.bam"
    filt_bam="${intermediate_dir}/${file_name}_mapq_30.bam"
    sample_counts="${output_dir}/${file_name}_miRNA_counts.txt"

    fastqc "$file" \
        -q \
        -o \
        "$intermediate_dir" &> /dev/null

    umi_tools extract \
        --extract-method=regex \
        --bc-pattern='.+(?P<discard_1>AACTGTAGGCACCATCAAT){s<=2}(?P<umi_1>.{12})(?P<discard_2>.+)' \
        -I "$file" \
        -S "$umi_extracted_file" \
        -L "${log_dir}/${file_name}_umi_extract.log"

    cutadapt \
        --minimum-length 18 \
        --poly-a \
        --quality-cutoff 20,20 \
        --output "$trimmed_file" \
        "$umi_extracted_file" > "${log_dir}/${file_name}_cutadapt.log"
    
    fastqc "$trimmed_file" \
        -q \
        -o \
        "$intermediate_dir" &> /dev/null
    
    awk 'NR%4==2 { print length($0) }' "$trimmed_file" \
        | sort -n | uniq -c | awk '{ print $2"\t"$1 }' \
        > "${output_dir}/stats/${file_name}_length_dist.tsv"

    bwa aln \
        -n 1 \
        -o 0 \
        -e 0 \
        -l 8 \
        -k 0 \
        -t 1 \
        -f "$bwa_output" \
        "$bwa_index" \
        "$trimmed_file" 2> "${log_dir}/${file_name}_bwa_aln.log"

    bwa samse \
        -f "$bwa_sam" \
        "$bwa_index" \
        "$bwa_output" \
        "$trimmed_file" 2> "${log_dir}/${file_name}_bwa_samse.log"

    samtools view -bS "$bwa_sam" | samtools sort -o "$bwa_bam"
    samtools index "$bwa_bam"

    samtools flagstat "$bwa_bam" > "${log_dir}/${file_name}_flagstat.txt"

    samtools view -b -q 30 "$bwa_bam" > "$filt_bam"
    samtools index "$filt_bam"

    umi_tools dedup \
        -I "$bwa_bam" \
        --log2stderr \
        -S "$dedup_bam" &> /dev/null

    samtools index "$dedup_bam"

    samtools idxstats "$dedup_bam" | awk 'NR > 1 { print prev } { prev = $1"\t"$3 }' > "$sample_counts"

    raw=$(( $(wc -l < "$file") / 4 ))
    umi=$(( $(wc -l < "$umi_extracted_file") / 4 ))
    trimmed=$(( $(wc -l < "$trimmed_file") / 4 ))
    mapped=$(samtools view -c -F 4 "$bwa_bam")
    passq=$(samtools view -c "$filt_bam")
    deduped=$(samtools view -c "$dedup_bam")
    in_range=$(awk 'NR%4==2 { l=length($0); if (l>=18 && l<=25) n++ } END { print n+0 }' "$trimmed_file")
    detected=$(awk '$2 >= 5' "$sample_counts" | wc -l | tr -d ' ')
 
    pct() { awk -v a="$1" -v b="$2" 'BEGIN { if (b==0) print "NA"; else printf "%.2f", 100*a/b }'; }
 
    printf '%s\t%d\t%d\t%d\t%s\t%d\t%s\t%d\t%s\t%d\n' \
        "$file_name" "$raw" "$umi" "$trimmed" \
        "$(pct "$in_range" "$trimmed")" \
        "$mapped" "$(pct "$mapped" "$trimmed")" \
        "$passq" "$(pct "$passq" "$trimmed")" \
        "$deduped" \
        > "${output_dir}/stats/${file_name}_stats.tsv"

    echo "Pipeline completed successfully for ${file_name}!"
}

export -f process_fastq

find "$fastq_dir" -name "*.fastq" -o -name "*.fq" | parallel  process_fastq {} "$output_dir" "$bwa_index"

{
  echo -e "Sample\tRaw reads\tUMI extracted reads\tTrimmed reads\t% reads within 18-25nt\tMapped reads\t% mapped\tMAPQ >= 30\t% MAPQ >= 30\tUnique molecules"
  cat "${output_dir}"/stats/*_stats.tsv
} > "${output_dir}/summary_stats.tsv"
 
echo "Summary written to ${output_dir}/summary_stats.tsv"