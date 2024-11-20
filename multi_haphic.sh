#!/bin/bash
# HapHiC: Reference-independent, allele-aware scaffolding tool
# Modified from https://github.com/zengxiaofei/HapHiC/

# Default thread count
threads=8

# Path to Juicer tools (can be overridden by setting JUICER_TOOL_PATH environment variable)
juicer_tool="/home/xuzhang/software/HapHiC/utils/juicer_tools.1.9.9_jcuda.0.8.jar"

# Usage function
usage() {
    echo "Usage: $0 [-t threads] [-i single_id | -l ids_file]"
    echo "  -t    Number of threads to use (default: 8)"
    echo "  -i    Single ID to process"
    echo "  -l    File containing a list of IDs"
    exit 1
}

# Check if required tools are available
for tool in bwa samtools filter_bam haphic seqkit juicer java; do
    if ! command -v $tool &>/dev/null; then
        echo "Error: $tool is not installed or not in PATH." >&2
        exit 1
    fi
done

# Parse command-line options
single_id=""
ids_file=""

while getopts "t:i:l:" opt; do
    case ${opt} in
        t ) threads=$OPTARG ;;
        i ) single_id=$OPTARG ;;
        l ) ids_file=$OPTARG ;;
        * ) usage ;;
    esac
done
shift $((OPTIND - 1))

# Validate input options
if [[ -z $single_id && -z $ids_file ]] || [[ -n $single_id && -n $ids_file ]]; then
    usage
fi

# Prepare IDs list
if [[ -n $single_id ]]; then
    ids=($single_id)
elif [[ -f $ids_file ]]; then
    mapfile -t ids < "$ids_file"
else
    echo "Error: File $ids_file does not exist."
    exit 1
fi

# Set working directory
wd=$(pwd)

# Process each ID
for id in "${ids[@]}"; do
    log_file="${wd}/${id}/${id}_haphic.log"

    # Validate working directory for ID
    id_dir="${wd}/${id}"
    if [[ ! -d $id_dir ]]; then
        echo "Error: Directory $id_dir does not exist. Skipping $id." >> "$log_file"
        continue
    fi

    pushd "$id_dir" > /dev/null || continue
	echo "$id - $(date)" &&  echo "$id - $(date)" | tee -a "$log_file"

    # Step 1: Index reference assemblies
    for hap in hap1 hap2; do
        ref_file="asm_hifiasm_hic/${id}.hifiasm.hic.${hap}.fa"
        bwt_file="${ref_file}.bwt"
        if [[ ! -f "$bwt_file" || ! -s "$bwt_file" ]]; then
            echo "Indexing reference for $hap..." >> "$log_file"
            bwa index "$ref_file" >> "$log_file" 2>&1 &
        else
            echo "Reference already indexed for $hap, skipping." >> "$log_file"
        fi
    done
    wait

    # Step 2: Align Hi-C reads and filter BAM files
    bam_status=$(samtools quickcheck ${id}.hap*.HiC.filtered.bam 2>/dev/null && echo "valid" || echo "invalid")
    if [[ $bam_status == "valid" ]]; then
        echo "Filtered BAM files are valid for $id. Skipping alignment." >> "$log_file"
    else
        for hap in hap1 hap2; do
            { 
                ref_file="asm_hifiasm_hic/${id}.hifiasm.hic.${hap}.fa"
                bam_file="${id}.${hap}.HiC.bam"
                filtered_bam="${id}.${hap}.HiC.filtered.bam"

                if [[ ! -f "$bam_file" || ! -s "$bam_file" ]]; then
                    echo "Aligning reads to generate BAM file for $hap..." >> "$log_file"
                    bwa mem -5SP -t "$threads" "$ref_file" hic/*_1.fastq.gz hic/*_2.fastq.gz | \
                    samblaster | samtools view -@ "$threads" -b -F 3340 -o "$bam_file" >> "$log_file" 2>&1
                else
                    echo "BAM file already exists for $hap, skipping alignment." >> "$log_file"
                fi

                if [[ ! -f "$filtered_bam" || ! -s "$filtered_bam" ]]; then
                    echo "Filtering BAM file for $hap..." >> "$log_file"
                    filter_bam "$bam_file" 1 --nm 3 --threads "$threads" | \
                    samtools view -b -@ "$threads" -o "$filtered_bam" >> "$log_file" 2>&1
                else
                    echo "Filtered BAM file already exists for $hap, skipping filtering." >> "$log_file"
                fi
                if [[ -s "$filtered_bam" ]]; then
                    rm "$bam_file" && echo "Removed HiC.bam file to save disk space" >> "$log_file" 2>&1
                fi
            } & 
        done
        wait
    fi

    # Step 3: Run HapHiC scaffolding
    nchr_file="nchr"
    if [[ ! -f $nchr_file ]]; then
        echo "Error: Chromosome count file ($nchr_file) is missing. Skipping HapHiC." >> "$log_file"
        popd > /dev/null
        continue
    fi

    for hap in hap1 hap2; do
        chr_dir="chr_$hap"
        ref_file="asm_hifiasm_hic/${id}.hifiasm.hic.${hap}.fa"
        gfa_files="${id}.hifiasm.hic.hap1.p_ctg.gfa,${id}.hifiasm.hic.hap2.p_ctg.gfa"
        if [[ ! -d "$chr_dir" || ! -s "$chr_dir/${id}_chr_$hap.fa" ]]; then
            rm -rf "$chr_dir"
            echo "Running HapHiC pipeline for $hap..." >> "$log_file"
            haphic pipeline "$ref_file" "${id}.${hap}.HiC.filtered.bam" $(cat $nchr_file) \
            --outdir "$chr_dir" --gfa "$gfa_files" --threads $threads >> "$log_file" 2>&1 &
        else
            echo "HapHiC pipeline already completed for $hap, skipping." >> "$log_file"
        fi
    done
    wait

    # Step 4: Extract scaffolding results
    for hap in hap1 hap2; do
        chr_dir="chr_$hap"
        if [[ -d "$chr_dir" ]]; then
            (
                cd "$chr_dir" || { echo "Failed to enter $chr_dir"; exit 1; }

                echo "Processing $hap started at $(date)" >> "$log_file"

                # Extract scaffolds
                if seqkit grep -r -p ^group 04.build/scaffolds.fa >"${id}_chr_$hap.fa"; then
                    echo "Scaffolds extracted for $hap" >> "$log_file"
                    seqkit stat -Ta "${id}_chr_$hap.fa" >"${id}_chr_$hap.stat"
                else
                    echo "Failed to extract scaffolds for $hap" >> "$log_file"
                    continue
                fi

                # Generate Juicer TXT
                if juicer pre -a -o ${id}_${hap} \
                    "$id_dir/${id}.${hap}.HiC.filtered.bam" \
                    04.build/scaffolds.agp \
                    "$id_dir/asm_hifiasm_hic/${id}.hifiasm.hic.${hap}.fa.fai" >> "${id}_${hap}_juicer.log" 2>&1; then
                    echo "Juicer TXT generation complete for $hap" >> "$log_file"
                else
                    echo "Juicer TXT generation failed for $hap" >> "$log_file"
                    continue
                fi

                # Generate Juicer HIC
                if java -Xmx30G -jar ${juicer_tool} \
                    pre ${id}_${hap}.txt ${id}_${hap}.hic.part \
                    <(grep PRE_C_SIZE ${id}_${hap}_juicer.log | awk '{print $2" "$3}') >> "$log_file" 2>&1; then
                    mv "${id}_${hap}.hic.part" "${id}_${hap}.hic"
                    echo "Generated .hic file for $hap: ${id}_${hap}.hic" >> "$log_file"
                else
                    echo "Failed to generate .hic file for $hap" >> "$log_file"
                fi
            ) &
        else
            echo "Skipping $hap as $chr_dir does not exist" >> "$log_file"
        fi
    done
    wait
    echo "Processing complete for $id." | tee -a "$log_file"
    popd > /dev/null
done

echo "All tasks completed at $(date)."
