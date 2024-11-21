#!/bin/bash

# Function to process a single sample ID
process_sample() {
  local i=$1
  local wd=$(pwd)
  local base_dir=$wd/${i}

  if [ ! -d "$base_dir" ]; then
      echo "Directory $base_dir does not exist. Skipping $i."
      return
  fi

  cd $base_dir || return

  echo "Preparing JBAT files for $i with juicer pre..."
  juicer pre -a -o chr_hap1/${i}_hap1 ${i}.hap1.HiC.filtered.bam chr_hap1/04.build/scaffolds.agp asm_hifiasm_hic/${i}.hifiasm.hic.hap1.fa.fai >${i}_hap1_juicer.log 2>&1 &
  juicer pre -a -o chr_hap2/${i}_hap2 ${i}.hap2.HiC.filtered.bam chr_hap2/04.build/scaffolds.agp asm_hifiasm_hic/${i}.hifiasm.hic.hap2.fa.fai >${i}_hap2_juicer.log 2>&1 &

  wait
  echo "Generating txt done for $i at $(date)"

  local juicer_tool="/home/xuzhang/software/HapHiC/utils/juicer_tools.1.9.9_jcuda.0.8.jar"

  echo "Generating .hic files for visualization..."
  cd $base_dir/chr_hap1
  (java -jar -Xmx32G $juicer_tool \
    pre ${i}_hap1.txt ${i}_hap1.hic.part \
    <(cat ../${i}_hap1_juicer.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) \
    && (mv ${i}_hap1.hic.part ${i}_hap1.hic) &

  cd $base_dir/chr_hap2
  (java -jar -Xmx32G $juicer_tool \
    pre ${i}_hap2.txt ${i}_hap2.hic.part \
    <(cat ../${i}_hap2_juicer.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) \
    && (mv ${i}_hap2.hic.part ${i}_hap2.hic) &

  wait
  cd $base_dir
  echo "Process complete! Hi-C files ready for visualization in Juicebox."

  echo "Cleaning up .txt files if .hic files exist..."
  if [ -f "chr_hap1/${i}_hap1.hic" ]; then
      rm -f "chr_hap1/${i}_hap1.txt"
      echo "Deleted chr_hap1/${i}_hap1.txt"
  else
      echo "chr_hap1/${i}_hap1.hic does not exist."
  fi

  if [ -f "chr_hap2/${i}_hap2.hic" ]; then
      rm -f "chr_hap2/${i}_hap2.txt"
      echo "Deleted chr_hap2/${i}_hap2.txt"
  else
      echo "chr_hap2/${i}_hap2.hic does not exist."
  fi

  # Return to the working directory
  cd $wd
}

# Main script logic
if [ "$1" == "-l" ]; then
  # Process a list of IDs
  if [ -z "$2" ]; then
      echo "Usage: $0 -l <id_list_file>"
      exit 1
  fi

  id_list_file=$2
  if [ ! -f "$id_list_file" ]; then
      echo "File $id_list_file does not exist. Exiting."
      exit 1
  fi

  while IFS= read -r id; do
      [ -z "$id" ] && continue  # Skip empty lines
      process_sample "$id"
  done < "$id_list_file"
else
  # Process a single ID
  if [ -z "$1" ]; then
      echo "Usage: $0 <sample_id> or $0 -l <id_list_file>"
      exit 1
  fi

  process_sample "$1"
fi
