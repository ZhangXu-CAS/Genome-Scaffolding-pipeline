# Genome-Scaffolding-pipeline

A pipeline used for scaffolding two haplotypes from hifiasm to chromosome-level with hic data.

  Options:
    -h, --help  show this help message and exit
    -sp SP      Sample prefix for output file naming
    -hap1 HAP1  Input GFA file for haplotype 1
    -hap2 HAP2  Input GFA file for haplotype 2
    -hic1 HIC1  Hi-C read 1 FASTQ file
    -hic2 HIC2  Hi-C read 2 FASTQ file
    -t T        Number of threads for multi-threaded tasks
    -nchr NCHR  Number of chromosomes for the HapHiC pipeline

##
  python haphic_for2haps.py  -sp sample_prefix -hap1 hap1.gfa -hap2 hap2.gfa -hic1 hic_R1.fastq.gz -hic2 hic_R2.fastq.gz -t threads -nchr num_chromosomes

  Outputs:
    sp_GFA_to_FASTA: Convert GFA to FASTA and generate stats
    sp_Filtered_BAM: Index and align Hi-C data to the assembly
    sp_HapHiC_Resul: HapHiC scaffolding results
