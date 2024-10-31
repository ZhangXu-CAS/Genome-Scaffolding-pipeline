# Genome-Scaffolding-pipeline

A HapHiC pipeline used for scaffolding two haplotypes from hifiasm to chromosomal level with Hi-C data.

## Requirements:
- **haphic**
- **gfatools**
- **seqkit**
- **samblaster**

## Options:

- **-h, --help**:  Show this help message and exit
- **-sp   SP**:  Sample prefix for output file naming
- **-hap1   HAP1**: Input GFA file for haplotype 1
- **-hap2   HAP2**: Input GFA file for haplotype 2
- **-hic1   HIC1**: Hi-C read 1 FASTQ file
- **-hic2   HIC2**: Hi-C read 2 FASTQ file
- **-t   T**: Number of threads for multi-threaded tasks
- **-nchr   NCHR**: Number of chromosomes for the HapHiC pipeline
- **-mq** type=int Mapping quality filter (default: 1)
- **-nm** type=int Edit distance filter (default: 3)

## Output:
-  sp_GFA_to_FASTA: Convert GFA to FASTA and generate stats
-  sp_Filtered_BAM: Index and align Hi-C data to the assembly
-  sp_HapHiC_Resul: HapHiC scaffolding results

  
## Example:

```bash
python haphic_for2haps.py -sp sample_prefix -hap1 hap1.gfa -hap2 hap2.gfa -hic1 hic_R1.fastq.gz -hic2 hic_R2.fastq.gz -t threads -nchr num_chromosomes

```
## References:

### **Haphic**
- Xiaofei Zeng, Zili Yi, Xingtan Zhang, Yuhui Du, Yu Li, Zhiqing Zhou, Sijie Chen, Huijie Zhao, Sai Yang, Yibin Wang, Guoan Chen. Chromosome-level scaffolding of haplotype-resolved assemblies using Hi-C data without reference genomes. Nature Plants, 10:1184-1200. doi: https://doi.org/10.1038/s41477-024-01755-3
