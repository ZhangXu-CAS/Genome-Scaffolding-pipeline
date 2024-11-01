import os
import sys
import subprocess
import argparse
import logging
from datetime import datetime

# Define argument parser
def parse_args():
    parser = argparse.ArgumentParser(description="This script scaffolds two haplotypes from hifiasm using Hi-C data.")
    parser.add_argument("-sp", required=True, help="Sample prefix for output file naming")
    parser.add_argument("-hap1", required=True, help="Input GFA file for haplotype 1")
    parser.add_argument("-hap2", required=True, help="Input GFA file for haplotype 2")
    parser.add_argument("-hic1", required=True, help="Hi-C read 1 FASTQ file")
    parser.add_argument("-hic2", required=True, help="Hi-C read 2 FASTQ file")
    parser.add_argument("-t", type=int, required=True, help="Number of threads for multi-threaded tasks")
    parser.add_argument("-nchr", type=int, required=True, help="Number of chromosomes for the HapHiC pipeline")
    parser.add_argument("-mq", type=int, default=1, help="Mapping quality filter (default: 1)")
    parser.add_argument("-nm", type=int, default=3, help="Edit distance filter (default: 3)")
    return parser.parse_args()

# Set up logging
def setup_logging(sp):
    log_file = f"{sp}.haphic.2hap.log"
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info("Log file created: %s", log_file)

# Define function for running shell commands with error checking
def run_command(command, description, output_file):
    if os.path.exists(output_file):
        logging.info("%s already completed.", description)
    else:
        logging.info("Running: %s", description)
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info("Command output: %s", result.stdout.decode())
        if result.returncode != 0:
            logging.error("Error in %s: %s", description, result.stderr.decode())
            sys.exit(1)

def main():
    args = parse_args()
    setup_logging(args.sp)

    # Set up output directories
    output_dirs = [f"{args.sp}_GFA_to_FASTA", f"{args.sp}_Filtered_BAM", f"{args.sp}_HapHiC_Result"]
    for directory in output_dirs:
        os.makedirs(directory, exist_ok=True)

    # Convert GFA to FASTA and generate stats
    logging.info("Converting GFA to FASTA and generating stats...")
    run_command(f"gfatools gfa2fa {args.hap1} > {args.sp}_GFA_to_FASTA/{args.sp}.hap1.fa", "GFA to FASTA for hap1", f"{args.sp}_GFA_to_FASTA/{args.sp}.hap1.fa")
    run_command(f"samtools faidx {args.sp}_GFA_to_FASTA/{args.sp}.hap1.fa", "samtools faidx for hap1", f"{args.sp}_GFA_to_FASTA/{args.sp}.hap1.fa.fai")
    run_command(f"gfatools gfa2fa {args.hap2} > {args.sp}_GFA_to_FASTA/{args.sp}.hap2.fa", "GFA to FASTA for hap2", f"{args.sp}_GFA_to_FASTA/{args.sp}.hap2.fa")
    run_command(f"samtools faidx {args.sp}_GFA_to_FASTA/{args.sp}.hap2.fa", "samtools faidx for hap2", f"{args.sp}_GFA_to_FASTA/{args.sp}.hap2.fa.fai")
    run_command(f"seqkit stat -Ta {args.sp}_GFA_to_FASTA/{args.sp}.hap1.fa {args.sp}_GFA_to_FASTA/{args.sp}.hap2.fa > {args.sp}_GFA_to_FASTA/{args.sp}.2hap.stat", "Seqkit stats", f"{args.sp}_GFA_to_FASTA/{args.sp}.2hap.stat")

    # Index and align Hi-C data to the assembly
    logging.info("Indexing FASTA files for BWA and aligning Hi-C data...")
    run_command(f"bwa index {args.sp}_GFA_to_FASTA/{args.sp}.hap1.fa", "BWA index for hap1", f"{args.sp}_GFA_to_FASTA/{args.sp}.hap1.fa.bwt")
    run_command(f"bwa index {args.sp}_GFA_to_FASTA/{args.sp}.hap2.fa", "BWA index for hap2", f"{args.sp}_GFA_to_FASTA/{args.sp}.hap2.fa.bwt")

    # Align Hi-C reads to hap1 and hap2 separately
    bwa_mem_command_hap1 = (
        f"bwa mem -5SP -t {args.t} {args.sp}_GFA_to_FASTA/{args.sp}.hap1.fa {args.hic1} {args.hic2} | "
        "samblaster | "
        f"samtools view -@ {args.t} -S -h -b -F 3340 -o {args.sp}_Filtered_BAM/{args.sp}.hap1.HiC.bam"
    )
    bwa_mem_command_hap2 = (
        f"bwa mem -5SP -t {args.t} {args.sp}_GFA_to_FASTA/{args.sp}.hap2.fa {args.hic1} {args.hic2} | "
        "samblaster | "
        f"samtools view -@ {args.t} -S -h -b -F 3340 -o {args.sp}_Filtered_BAM/{args.sp}.hap2.HiC.bam"
    )
    run_command(bwa_mem_command_hap1, "Align Hi-C to hap1", f"{args.sp}_Filtered_BAM/{args.sp}.hap1.HiC.bam")
    run_command(bwa_mem_command_hap2, "Align Hi-C to hap2", f"{args.sp}_Filtered_BAM/{args.sp}.hap2.HiC.bam")

    # Filter alignments
    logging.info("Filtering alignments for MAPQ and NM...")
    filter_command_hap1 = (
        f"filter_bam {args.sp}_Filtered_BAM/{args.sp}.hap1.HiC.bam {args.mq} --nm {args.nm} --threads {args.t} | "
        f"samtools view -@ {args.t} -b -o {args.sp}_Filtered_BAM/{args.sp}.hap1.HiC.filtered.bam"
    )
    filter_command_hap2 = (
        f"filter_bam {args.sp}_Filtered_BAM/{args.sp}.hap2.HiC.bam {args.mq} --nm {args.nm} --threads {args.t} | "
        f"samtools view -@ {args.t} -b -o {args.sp}_Filtered_BAM/{args.sp}.hap2.HiC.filtered.bam"
    )
    run_command(filter_command_hap1, "Filter alignments for hap1", f"{args.sp}_Filtered_BAM/{args.sp}.hap1.HiC.filtered.bam")
    run_command(filter_command_hap2, "Filter alignments for hap2", f"{args.sp}_Filtered_BAM/{args.sp}.hap2.HiC.filtered.bam")
    os.remove(f"{args.sp}_Filtered_BAM/{args.sp}.hap1.HiC.bam")
    os.remove(f"{args.sp}_Filtered_BAM/{args.sp}.hap2.HiC.bam")

    # Run HapHiC scaffolding pipeline
    logging.info("Running HapHiC scaffolding pipeline...")
    haphic_command = (
        f"haphic pipeline {args.sp}_GFA_to_FASTA/{args.sp}.hap{{}}.fa {args.sp}_Filtered_BAM/{args.sp}.hap{{}}.HiC.filtered.bam {args.nchr} "
        f"--outdir {args.sp}_HapHiC_Result/chr_hap{{}} --gfa {args.hap1},{args.hap2} --threads {args.t}"
    )
    run_command(haphic_command.format(1), "HapHiC pipeline for hap1", f"{args.sp}_HapHiC_Result/chr_hap1")
    run_command(haphic_command.format(2), "HapHiC pipeline for hap2", f"{args.sp}_HapHiC_Result/chr_hap2")

    # Process HapHiC results
    for hap in ["hap1", "hap2"]:
        logging.info("Processing %s results...", hap)
        scaffold_fa = f"{args.sp}_HapHiC_Result/{args.sp}_chr_{hap}.fa"
        scaffolds_dir = f"{args.sp}_HapHiC_Result/chr_{hap}/04.build"
        run_command(f"seqkit grep -r -p ^group {scaffolds_dir}/scaffolds.fa > {scaffold_fa}", f"Scaffold for {hap}", scaffold_fa)
    
    logging.info("All processes completed successfully!")

if __name__ == "__main__":
    main()