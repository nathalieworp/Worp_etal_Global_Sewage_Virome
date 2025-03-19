# Snakemake pipeline for analysis of global urban sewage metagenomic samples.
# The pipeline is designed to run sequentially per batch and method. 
#
# Required folder structure (must be created manually before running the pipeline):
#
# /mnt/viro0002/workgroups_projects/Global_Sewage/batches/{batch}/{method}/
# ├── workflow_metaspades/  # Output files from the pipeline will be stored here
# ├── raw_data/             # Raw sequencing reads should be placed here
#
# Example:
# /mnt/viro0002/workgroups_projects/Global_Sewage/batches/batch1/capture/
# ├── workflow_metaspades/
# ├── raw_data/  (contains {sample}_R1_001.fastq.gz and {sample}_R2_001.fastq.gz)
#
# **INSTRUCTIONS:**
# 1. Manually create the required folders** (`workflow_metaspades/` and `raw_data/`).
# 2. Ensure raw sequencing reads are in `raw_data/` inside each batch/method folder.
# 3. Run the pipeline from inside `workflow_metaspades/`
# 4. Output files will be generated inside `workflow_metaspades/`.
#
# **Note:** Running the pipeline outside `workflow_metaspades/` may cause errors.
import os

FILES, = glob_wildcards("raw_data/{sample}_R1_001.fastq.gz")
batches = ["batch1", "batch2", "batch3","batch4","batch5","batch6","batch7"]
methods = ["metagenomics","capture"]

rule all:
    input:
        expand("{sample}/assembly/contigs.fasta",sample=FILES),
        expand("{sample}/{sample}_bam_stat.tsv",sample=FILES),
        expand("{sample}/{sample}_annotation.tsv",sample=FILES),
        expand("{sample}/completed_{sample}_annotation.tsv", sample=FILES)

rule unzip_file:
    input:
        R1 = "raw_data/{sample}_R1_001.fastq.gz",
        R2 = "raw_data/{sample}_R2_001.fastq.gz"
    output:
        R1 = temp("{sample}/{sample}_R1_001.fastq"),
        R2 = temp("{sample}/{sample}_R2_001.fastq")
    threads: 4
    shell:
        """
        cat {input.R1} | unpigz -k -p {threads} > {output.R1}
        cat {input.R2} | unpigz -k -p {threads} > {output.R2}
        """
#count the number of reads in sample expr $(cat file.fastq | wc -l) / 4
#CD-HIT-dup identifies duplicates from single or paired Illumina reads
#first deduplicate before QC (else you can have trimmed sequences of different length)

rule dedup_raw:
    input:
        R1 = "{sample}/{sample}_R1_001.fastq",
        R2 = "{sample}/{sample}_R2_001.fastq"
    output:
        R1 = "{sample}/dedup/{sample}_R1_001.fastq.gz",
        R2 = "{sample}/dedup/{sample}_R2_001.fastq.gz"
    threads: 4
    shell:
        """
        cd-hit-dup -u 150 -i {input.R1} -i2 {input.R2} \
        -o {wildcards.sample}/dedup/R1.tmp -o2 {wildcards.sample}/dedup/R2.tmp
        pigz -p {threads} {wildcards.sample}/dedup/R1.tmp {wildcards.sample}/dedup/R2.tmp
        mv {wildcards.sample}/dedup/R1.tmp.gz {output.R1}
        mv {wildcards.sample}/dedup/R2.tmp.gz {output.R2}
        rm {wildcards.sample}/dedup/*.clstr
        """

rule QC_after_dedup:
    input:
        R1 = "{sample}/dedup/{sample}_R1_001.fastq.gz",
        R2 = "{sample}/dedup/{sample}_R2_001.fastq.gz"
    output:
        R1 = "{sample}/dedup_qc/{sample}_R1_001.fastq.gz",
        R2 = "{sample}/dedup_qc/{sample}_R2_001.fastq.gz",
        S = "{sample}/dedup_qc/{sample}_S_001.fastq.gz",
        failed = "{sample}/dedup_qc/{sample}_fail_001.fastq.gz"
    threads: 1
    shell:
        """
        fastp -i {input.R1} -I {input.R2} \
        -o {output.R1} -O {output.R2} \
        --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        --unpaired1 {output.S} --unpaired2 {output.S} --failed_out {output.failed} \
        -l 30 -Q -y \
        --cut_right \
        --cut_right_window_size 5 \
        --cut_right_mean_quality 25 \
        -w {threads} \
        -j {wildcards.sample}/dedup_qc/qc_report.json -h {wildcards.sample}/dedup_qc/qc_report.html
        """
        
#for PE input, if read1 passed QC but read2 not, it will be written to unpaired1. Default is to discard it.
# file S contains singleton reads. If --unpaired2 is same as --unpaired1 (default mode), both unpaired reads will be written to this same file. (string [=])

rule filter_human:
    input:
        R1 = "{sample}/dedup_qc/{sample}_R1_001.fastq.gz",
        R2 = "{sample}/dedup_qc/{sample}_R2_001.fastq.gz",
        S = "{sample}/dedup_qc/{sample}_S_001.fastq.gz"
    output:
        R1 = "{sample}/filtered/{sample}_R1_001.fastq.gz",
        R2 = "{sample}/filtered/{sample}_R2_001.fastq.gz",
	S = "{sample}/filtered/{sample}_S_001.fastq.gz"
    threads: 4
    shell:
        """
        bwa mem -aY -t {threads} /mnt/viro0002/home/r059576/Sewage_study_nath/genome_GCF/GCF_000001405.26_GRCh38_genomic.fna {input.R1} {input.R2} | \
        samtools fastq -c 6 -f 4 -1 {output.R1} -2 {output.R2} -s /dev/null - # - standardin --> waar de standardin erin moet; singlets worden weggegooid (waarvan 1 read mapt tegen humaan genoom) 
        bwa mem -aY -t {threads} /mnt/viro0002/home/r059576/Sewage_study_nath/genome_GCF/GCF_000001405.26_GRCh38_genomic.fna {input.S} | \
        samtools fastq -c 6 -f 4 - | gzip > {output.S}  
        """

# -aY . Y is use softclipping -a is output all alignments
# the output of bwa mem is a sam file that is piped directly into samtools in order to not save the memory consuming sam files. The samtools fastq option -c representes the level of compression when writing .gz fastq files. The output bam/sam file contains information about mapped and unmapped reads. 
# -f 4 : then you obtain the unmapped (which are not human reads) reads and keep these. 
#-s option: If a singleton file is specified using the -s option then only paired sequences will be output for categories 1 and 2; 1 and 2. This can be used to prepare fastq files for programs that cannot handle a mixture of paired and singleton reads.

rule assemble_filtered:
    input:
        R1 = "{sample}/filtered/{sample}_R1_001.fastq.gz",
        R2 = "{sample}/filtered/{sample}_R2_001.fastq.gz",
        S = "{sample}/filtered/{sample}_S_001.fastq.gz"
    output:
        "{sample}/assembly/contigs.fasta"
    params:
        spades_folder = "{sample}/assembly"
    threads: 16
    shell:
        """
        spades.py -t {threads} \
        --meta \
        -o {params.spades_folder} \
        -1 {input.R1} \
        -2 {input.R2} \
        -s {input.S}
        """
#de novo assembly of contigs 
#file with forward paired-end reads + file with reverse paired end reads
# -s file with unpaired reads 

rule blastx_assembled:
    input:
        ancient("{sample}/assembly/contigs.fasta")
    output:
        "{sample}/{sample}_annotation.tsv"
    threads: 48
    shell:
        """
        diamond blastx \
        -q {input} \
	-d /mnt/viro0002/workgroups_projects/Bioinformatics/DB/DIAMOND_db_2022-11-23/2022-11-23_nr_db.dmnd \
        -o {output} \
        -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids \
        --taxonmap /mnt/viro0002/workgroups_projects/Bioinformatics/DB/taxdump+prot2accetaxid_2022-11-23/prot.accession2taxid.FULL.gz \
        --threads {threads} \
        -b10.0
        """
rule map_reads_to_contigs:
    input:
        R1 = "{sample}/filtered/{sample}_R1_001.fastq.gz",
        R2 = "{sample}/filtered/{sample}_R2_001.fastq.gz",
        S = "{sample}/filtered/{sample}_S_001.fastq.gz",
        contigs = ancient("{sample}/assembly/contigs.fasta")
    output:
        "{sample}/{sample}_mappings.bam"
    threads: 4
    shell:
        """
        bwa index {input.contigs}
        bwa mem -Y -t {threads} {input.contigs} {input.R1} {input.R2} | samtools sort - > {wildcards.sample}/tmp_paired.bam
        bwa mem -Y -t {threads} {input.contigs} {input.S} | samtools sort - > {wildcards.sample}/tmp_singlets.bam
        samtools merge {output} {wildcards.sample}/tmp_paired.bam {wildcards.sample}/tmp_singlets.bam
        rm {wildcards.sample}/tmp_paired.bam {wildcards.sample}/tmp_singlets.bam
        """
rule create_coverage_file:
    input:
        "{sample}/{sample}_mappings.bam"
    output:
        "{sample}/{sample}_coverage.txt"
    threads: 1
    shell:
        """
        samtools view -bF2052 {input} | samtools depth -a -d 0 - > {output}
        """
        
# -b: output in de bam format, -F2025 Do no2 output alignments with any bits set in INT present in the FLAG field. INT = 2025 and means that alignments that are not mapped are not outputted. 
#samtools depth: Computes the depth at each position or region. -a Output all positions (including those with zero depth) -d 
# -d 0:  read at most 0 (INT) reads per input file. This means figures greater than INT may be reported in the output.Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. Passing zero for this option sets it to the highest possible value, effectively removing the depth limit.

rule create_readmap_file:
    input:
        "{sample}/{sample}_mappings.bam"
    output:
        "{sample}/{sample}_readmap.txt"
    threads: 1
    shell:
        """
        samtools view -F2052 {input} | cut -f1,3 > {output}
        """
#, do not output unmapped alignments and print the first and third field to output 

rule stat_mapped_files:
    input:
        "{sample}/{sample}_mappings.bam"

    output:
        "{sample}/{sample}_bam_stat.tsv"

    threads: 1
    shell:
        """
        samtools view -c -F2052 {input} | cut -f1,3 > {output}
        """

rule merge_results:
    input:
        annotation = "{sample}/{sample}_annotation.tsv",
        coverage = "{sample}/{sample}_coverage.txt",
        readmap = "{sample}/{sample}_readmap.txt"
    output:
        "{sample}/completed_{sample}_annotation.tsv"
    threads: 1
    script:
        "/mnt/viro0002/workgroups_projects/Global_Sewage/scripts/merge_results.R"


rule raw_stats:
    input:
        r1=expand("../raw_data/{sample}_R1_001.fastq.gz", sample=FILES),
        r2=expand("../raw_data/{sample}_R2_001.fastq.gz", sample=FILES)
    output:
        "../readstats_raw.tsv"
    shell:
        """
        seqkit stats {input.r1} {input.r2} | awk 'FNR==1 && NR!=1{{next}}{{print}}' > {output}
        """

rule dedup_stats:
    input:
        r1=expand("{sample}/dedup/{sample}_R1_001.fastq.gz", sample=FILES),
        r2=expand("{sample}/dedup/{sample}_R2_001.fastq.gz", sample=FILES)
    output:
        "../readstats_dedup.tsv"
    shell:
        """
        seqkit stats {input.r1} {input.r2} | awk 'FNR==1 && NR!=1{{next}}{{print}}' > {output}
        """
rule dedup_qc_stats:
    input:
        r1=expand("{sample}/dedup_qc/{sample}_R1_001.fastq.gz", sample=FILES),
        r2=expand("{sample}/dedup_qc/{sample}_R2_001.fastq.gz", sample=FILES)
    output:
        "../readstats_dedup_qc.tsv"
    shell:
        """
        seqkit stats {input.r1} {input.r2} | awk 'FNR==1 && NR!=1{{next}}{{print}}' > {output}
        """

rule filtered_stats:
    input:
        r1=expand("{sample}/filtered/{sample}_R1_001.fastq.gz", sample=FILES),
        r2=expand("{sample}/filtered/{sample}_R2_001.fastq.gz", sample=FILES)
    output:
        "../readstats_filtered.tsv"
    shell:
        """
        seqkit stats {input.r1} {input.r2} | awk 'FNR==1 && NR!=1{{next}}{{print}}' > {output}
        """

# This Snakemake pipeline was originally developed by David Nieuwenhuijse
# and modified/extended by Nathalie Worp
