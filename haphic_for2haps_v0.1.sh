#!/bin/bash
# This script is used to scaffold two haplotypes from hifiasm using Hi-C data.
# Draft by Xu Zhang and ChatGPT on Oct. 31, 2024

show_help() {
	echo ""
	echo "This script is used to scaffold two haplotypes from hifiasm using Hi-C data. "
    echo "Version: v.0.1"
    echo ""
    echo "Usage: bash haphic_for2haps_xx.sh -sp sample_prefix -hap1 hap1.gfa -hap2 hap2.gfa -hic1 hic_R1.fastq.gz -hic2 hic_R2.fastq.gz -t threads -nchr num_chromosomes [-mq map_quality] [-nm edit_distance]"
    echo ""
    echo "Options:"
    echo "  -sp    Sample prefix for output file naming (required)"
    echo "  -hap1  Input GFA file for haplotype 1 (required)"
    echo "  -hap2  Input GFA file for haplotype 2 (required)"
    echo "  -hic1  Hi-C read 1 FASTQ file (required)"
    echo "  -hic2  Hi-C read 2 FASTQ file (required)"
    echo "  -t     Number of threads to use for multi-threaded tasks (required)"
    echo "  -nchr  Number of chromosomes for the HapHiC pipeline (required)"
    echo "  -mq    Mapping quality filter for filter_bam (default: 1)"
    echo "  -nm    Edit distance filter for filter_bam (default: 3)"
    echo "  -v     Display version information"
    echo "  -h     Display this help message"
    exit 0
}

# Set default values for MQ and NM filters
mapq=1
nm=3

# Define version
version="v.0.1"

# Parse input arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -sp) sp="$2"; shift ;;
        -hap1) hap1="$2"; shift ;;
        -hap2) hap2="$2"; shift ;;
        -hic1) hic1="$2"; shift ;;
        -hic2) hic2="$2"; shift ;;
        -t) threads="$2"; shift ;;
        -nchr) nchr="$2"; shift ;;
        -mq) mapq="$2"; shift ;;
        -nm) nm="$2"; shift ;;
        -v) echo "Version: $version"; exit 0 ;;
        -h) show_help ;;
        *) echo "Unknown parameter passed: $1"; show_help ;;
    esac
    shift
done

# Ensure required parameters are provided
if [[ -z "$sp" || -z "$hap1" || -z "$hap2" || -z "$hic1" || -z "$hic2" || -z "$threads" || -z "$nchr" ]]; then
    echo "Error: Missing parameters."
    show_help
    exit 1
fi

# Define log file
log_file="${sp}_haphic_for2haps.log"

# Create output directories
mkdir -p "${sp}_GFA_to_FASTA" "${sp}_Filtered_BAM" "${sp}_HapHiC_Result"

# Define function for error checking
check_success() {
  if [ $? -ne 0 ]; then
    echo "Error in $1" | tee -a "$log_file"
    exit 1
  fi
}

# Convert GFA to FASTA and generate stats
echo "Converting GFA to FASTA and indexing..." | tee -a "$log_file"
gfatools gfa2fa "$hap1" > "${sp}_GFA_to_FASTA/${sp}.hap1.fa"
check_success "gfa2fa for hap1"
samtools faidx "${sp}_GFA_to_FASTA/${sp}.hap1.fa"
check_success "samtools faidx for hap1"

gfatools gfa2fa "$hap2" > "${sp}_GFA_to_FASTA/${sp}.hap2.fa"
check_success "gfa2fa for hap2"
samtools faidx "${sp}_GFA_to_FASTA/${sp}.hap2.fa"
check_success "samtools faidx for hap2"

seqkit stat -Ta "${sp}_GFA_to_FASTA/${sp}.hap1.fa" "${sp}_GFA_to_FASTA/${sp}.hap2.fa" > "${sp}_GFA_to_FASTA/${sp}.2hap.stat"
check_success "seqkit stat"

# Index and align Hi-C data to the assembly
echo "Indexing FASTA files for BWA..." | tee -a "$log_file"
bwa index "${sp}_GFA_to_FASTA/${sp}.hap1.fa" 
bwa index "${sp}_GFA_to_FASTA/${sp}.hap2.fa" 

check_success "bwa index"

echo "Aligning Hi-C data to the assemblies..." | tee -a "$log_file"
bwa mem -5SP -t "$threads" "${sp}_GFA_to_FASTA/${sp}.hap1.fa" "$hic1" "$hic2" | \
    samblaster | \
    samtools view -@ "$threads" -S -h -b -F 3340 -o "${sp}_Filtered_BAM/${sp}.hap1.HiC.bam" &
bwa mem -5SP -t "$threads" "${sp}_GFA_to_FASTA/${sp}.hap2.fa" "$hic1" "$hic2" | \
    samblaster | \
    samtools view -@ "$threads" -S -h -b -F 3340 -o "${sp}_Filtered_BAM/${sp}.hap2.HiC.bam" &
wait
check_success "bwa mem alignments"

# Filter alignments with MAPQ and NM
echo "Filtering alignments with MAPQ: $mapq and NM: $nm" | tee -a "$log_file"
filter_bam "${sp}_Filtered_BAM/${sp}.hap1.HiC.bam" "$mapq" --nm "$nm" --threads "$threads" | \
    samtools view -@ "$threads" -b -o "${sp}_Filtered_BAM/${sp}.hap1.HiC.filtered.bam" &
filter_bam "${sp}_Filtered_BAM/${sp}.hap2.HiC.bam" "$mapq" --nm "$nm" --threads "$threads" | \
    samtools view -@ "$threads" -b -o "${sp}_Filtered_BAM/${sp}.hap2.HiC.filtered.bam" &
wait
check_success "filter_bam"
rm "${sp}_Filtered_BAM/${sp}.hap1.HiC.bam" "${sp}_Filtered_BAM/${sp}.hap2.HiC.bam"

# Run HapHiC scaffolding pipeline
echo "Running HapHiC scaffolding pipeline..." | tee -a "$log_file"
haphic pipeline "${sp}_GFA_to_FASTA/${sp}.hap1.fa" "${sp}_Filtered_BAM/${sp}.hap1.HiC.filtered.bam" "$nchr" \
    --outdir "${sp}_HapHiC_Result/chr_hap1" --gfa "$hap1,$hap2" --threads "$threads" 2>>"$log_file" &
haphic pipeline "${sp}_GFA_to_FASTA/${sp}.hap2.fa" "${sp}_Filtered_BAM/${sp}.hap2.HiC.filtered.bam" "$nchr" \
    --outdir "${sp}_HapHiC_Result/chr_hap2" --gfa "$hap1,$hap2" --threads "$threads" 2>>"$log_file" &
wait
check_success "haphic pipeline"

# Process HapHiC scaffolding results
for hap in hap1 hap2; do
  echo "Processing $hap results..." | tee -a "$log_file"
  cd "${sp}_HapHiC_Result/chr_$hap" || exit
  seqkit grep -r -p ^group 04.build/scaffolds.fa > "${sp}_HapHiC_Result/${sp}_chr_$hap.fa"
  seqkit stat -Ta "${sp}_HapHiC_Result/${sp}_chr_$hap.fa" > "${sp}_HapHiC_Result/${sp}_chr_$hap.stat"
  haphic plot 04.build/scaffolds.raw.agp "../${sp}_Filtered_BAM/${sp}.$hap.HiC.filtered.bam" | tee -a "$log_file"
  cd - || exit
done

echo "Pipeline completed at: $(date)" | tee -a "$log_file"
