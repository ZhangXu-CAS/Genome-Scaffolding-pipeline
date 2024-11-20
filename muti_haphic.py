import os
import subprocess
import argparse
import logging
from multiprocessing import Pool

# Default thread count
threads = 8

# Path to Juicer tools (can be overridden by setting JUICER_TOOL_PATH environment variable)
juicer_tool = os.getenv("JUICER_TOOL_PATH", "/path/to/juicer_tools.1.9.9_jcuda.0.8.jar")

# Set up logging
logging.basicConfig(level=logging.INFO)

# Function to run a command and capture its output
def run_command(command, log_file):
    try:
        result = subprocess.run(command, shell=True, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logging.info(f"Command succeeded: {command}")
        with open(log_file, "a") as f:
            f.write(result.stdout)
        return result.stdout
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {command}")
        with open(log_file, "a") as f:
            f.write(f"Error: {e.stderr}\n")
        return None

# Validate if the required tools are installed
def check_tools():
    required_tools = ["bwa", "samtools", "filter_bam", "haphic", "seqkit", "juicer", "java"]
    for tool in required_tools:
        if shutil.which(tool) is None:
            raise EnvironmentError(f"Error: {tool} is not installed or not in PATH.")
    
# Function to handle ID processing
def process_id(id, threads, juicer_tool, wd):
    log_file = os.path.join(wd, id, f"{id}_haphic.log")
    id_dir = os.path.join(wd, id)

    if not os.path.isdir(id_dir):
        logging.error(f"Directory {id_dir} does not exist. Skipping {id}.")
        return

    os.chdir(id_dir)
    logging.info(f"Processing {id} - {os.getcwd()} at {id_dir}")

    # Step 1: Index reference assemblies
    for hap in ["hap1", "hap2"]:
        ref_file = f"asm_hifiasm_hic/{id}.hifiasm.hic.{hap}.fa"
        bwt_file = f"{ref_file}.bwt"
        if not os.path.isfile(bwt_file) or os.path.getsize(bwt_file) == 0:
            logging.info(f"Indexing reference for {hap}...")
            run_command(f"bwa index {ref_file}", log_file)

    # Step 2: Align Hi-C reads and filter BAM files
    bam_status = "valid" if subprocess.run(f"samtools quickcheck {id}.hap*.HiC.filtered.bam", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).returncode == 0 else "invalid"
    
    if bam_status == "valid":
        logging.info(f"Filtered BAM files are valid for {id}. Skipping alignment.")
    else:
        for hap in ["hap1", "hap2"]:
            ref_file = f"asm_hifiasm_hic/{id}.hifiasm.hic.{hap}.fa"
            bam_file = f"{id}.{hap}.HiC.bam"
            filtered_bam = f"{id}.{hap}.HiC.filtered.bam"
            
            if not os.path.isfile(bam_file) or os.path.getsize(bam_file) == 0:
                logging.info(f"Aligning reads to generate BAM file for {hap}...")
                bwa_command = f"bwa mem -5SP -t {threads} {ref_file} hic/*_1.fastq.gz hic/*_2.fastq.gz | samblaster | samtools view -@ {threads} -b -F 3340 -o {bam_file}"
                run_command(bwa_command, log_file)
            
            if not os.path.isfile(filtered_bam) or os.path.getsize(filtered_bam) == 0:
                logging.info(f"Filtering BAM file for {hap}...")
                filter_bam_command = f"filter_bam {bam_file} 1 --nm 3 --threads {threads} | samtools view -b -@ {threads} -o {filtered_bam}"
                run_command(filter_bam_command, log_file)
            
            if os.path.isfile(filtered_bam) and os.path.getsize(filtered_bam) > 0:
                os.remove(bam_file)
                logging.info(f"Removed HiC.bam file to save disk space for {hap}")

    # Step 3: Run HapHiC scaffolding
    nchr_file = "nchr"
    if not os.path.isfile(nchr_file):
        logging.error(f"Chromosome count file ({nchr_file}) is missing. Skipping HapHiC.")
        return

    for hap in ["hap1", "hap2"]:
        chr_dir = f"chr_{hap}"
        ref_file = f"asm_hifiasm_hic/{id}.hifiasm.hic.{hap}.fa"
        gfa_files = f"{id}.hifiasm.hic.hap1.p_ctg.gfa,{id}.hifiasm.hic.hap2.p_ctg.gfa"
        
        if not os.path.isdir(chr_dir) or not os.path.isfile(f"{chr_dir}/{id}_chr_{hap}.fa"):
            os.makedirs(chr_dir, exist_ok=True)
            logging.info(f"Running HapHiC pipeline for {hap}...")
            haphic_command = f"haphic pipeline {ref_file} {id}.{hap}.HiC.filtered.bam $(cat {nchr_file}) --outdir {chr_dir} --gfa {gfa_files} --threads {threads}"
            run_command(haphic_command, log_file)

    # Step 4: Extract scaffolding results
    for hap in ["hap1", "hap2"]:
        chr_dir = f"chr_{hap}"
        if os.path.isdir(chr_dir):
            scaffolds_file = f"{id}_chr_{hap}.fa"
            scaffolds_stat_file = f"{id}_chr_{hap}.stat"
            logging.info(f"Processing {hap}...")
            scaffolds_command = f"seqkit grep -r -p ^group 04.build/scaffolds.fa > {scaffolds_file}"
            run_command(scaffolds_command, log_file)
            seqkit_stat_command = f"seqkit stat -Ta {scaffolds_file} > {scaffolds_stat_file}"
            run_command(seqkit_stat_command, log_file)

            # Generate Juicer TXT
            juicer_txt_command = f"juicer pre -a -o {id}_{hap} {id_dir}/{id}.{hap}.HiC.filtered.bam 04.build/scaffolds.agp {id_dir}/asm_hifiasm_hic/{id}.hifiasm.hic.{hap}.fa.fai"
            run_command(juicer_txt_command, log_file)

            # Generate Juicer HIC
            java_command = f"java -Xmx30G -jar {juicer_tool} pre {id}_{hap}.txt {id}_{hap}.hic.part <(grep PRE_C_SIZE {id}_{hap}_juicer.log | awk '{{print $2\" \"$3}}')"
            run_command(java_command, log_file)

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="HapHiC: Reference-independent, allele-aware scaffolding tool")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads to use (default: 8)")
    parser.add_argument("-i", "--single_id", help="Single ID to process")
    parser.add_argument("-l", "--ids_file", help="File containing a list of IDs")
    args = parser.parse_args()

    # Validate input
    if not args.single_id and not args.ids_file:
        logging.error("Please specify either a single ID (-i) or a list of IDs (-l).")
        exit(1)

    if args.single_id and args.ids_file:
        logging.error("Please specify only one of -i or -l, not both.")
        exit(1)

    # Check if required tools are installed
    check_tools()

    # Prepare list of IDs
    if args.single_id:
        ids = [args.single_id]
    elif args.ids_file:
        with open(args.ids_file, "r") as f:
            ids = [line.strip() for line in f]

    # Get current working directory
    wd = os.getcwd()

    # Process each ID
    with Pool(threads) as pool:
        pool.starmap(process_id, [(id, args.threads, juicer_tool, wd) for id in ids])

    logging.info("All tasks completed.")