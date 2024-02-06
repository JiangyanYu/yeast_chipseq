"""
Author: Jiangyan Yu, modified from jsschrepping/bioinfo-base-image:jss_v0.0.3
Affiliation: LIMES, Uni Bonn
Aim: A simple Snakemake workflow to process paired-end stranded Chip-Seq. (modified from the ATAC-seq pipeline)
Date: 20240103
Run: snakemake --use-conda
Docker: docker run -it --rm --name Chip_preprocessing -v /path/to/data/:/data/ -v /path/to/index/:/data/bowtie2-index/genomefiles/ jsschrepping/bioinfo-base-image:jss_v0.0.3 /bin/bash
# bowtie2 index basename is the directory plus file name before .1.bt2

"""
#SAMPLES=["NG-31750_Zuo1MYC_2_lib648658_10157_2"]
SAMPLES=["NG-31750_Zuo1MYC_1_lib648657_10157_2","NG-31750_Zuo1MYC_2_lib648658_10157_2","NG-31750_Zuo1MYC_UV_1_lib648659_10157_2","NG-31750_Zuo1MYC_UV_2_lib648660_10157_2","NG-31750_Input_Zuo1MYC_lib648656_10157_2"]

READS=["1", "2"] # it is crucial to recognize sample and reads information for snakemake wildcard#

### Target rule ###

rule all:
    input:
        expand("output/qc/fastqc/{sample}_R{read}_fastqc.html", sample=SAMPLES, read=READS),
        expand("output/trimmed/{sample}_R{read}.fastq", sample=SAMPLES, read=READS),
        expand("output/qc/fastqc/{sample}_R{read}_trimmed_fastqc.html", sample=SAMPLES, read=READS),
        expand("output/mapped/{sample}.final.bam",sample=SAMPLES),
        expand("output/mapped/{sample}.final.bai",sample=SAMPLES),
        expand("output/mapped/{sample}.final.bw",sample=SAMPLES),
        expand("output/mapped/{sample}.final.flagstat",sample=SAMPLES),
        expand("output/mapped/{sample}.subsample.bam",sample=SAMPLES),
        expand("output/mapped/{sample}.subsample.bai",sample=SAMPLES),
        expand("output/mapped/{sample}.subsample.bw",sample=SAMPLES),
        expand("output/peaks/{sample}_peaks.xls",sample=SAMPLES),
        expand("output/peaks/{sample}_subsample_peaks.xls",sample=SAMPLES),
        "output/qc/multiqc_raw.html",
        "output/qc/multiqc.html"

### raw fastq ###

rule fastqc_raw:
    input:
        "/data/fastq/{sample}_R{read}.fastq"
    output:
        html="output/qc/fastqc/{sample}_R{read}_fastqc.html",
        zip="output/qc/fastqc/{sample}_R{read}_fastqc.zip"
    params: ""
    log:
        "output/logs/fastqc/{sample}_R{read}.log"
    wrapper:
        "0.34.0/bio/fastqc"

rule multiqc_raw:
    input:
        expand("output/qc/fastqc/{sample}_R{read}_fastqc.zip",zip,sample=SAMPLES,read=READS)
    output:
        "output/qc/multiqc_raw.html"
    log:
        "output/logs/multiqc_raw.log"
    params:
        ""
    wrapper:        
        "v3.3.3/bio/multiqc"

### Adapter trimming ###

rule trimmomatic_pe:
    input:
        r1="/data/fastq/{sample}_R1.fastq",
        r2="/data/fastq/{sample}_R2.fastq"
    output:
        r1="output/trimmed/{sample}_R1.fastq",
        r2="output/trimmed/{sample}_R2.fastq",
        # reads where trimming entirely removed the mate
        r1_unpaired="output/trimmed/{sample}_R1.unpaired.fastq",
        r2_unpaired="output/trimmed/{sample}_R2.unpaired.fastq"
    log:
        "output/logs/trimmomatic/{sample}.log"
    params:
        trimmer=["ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36"],
        extra="",
        compression_level="-9"
    threads:
        28
    wrapper:
        "0.35.0/bio/trimmomatic/pe"

rule fastqc_trimmed:
    input:
        "output/trimmed/{sample}_R{read}.fastq"
    output:
        html="output/qc/fastqc/{sample}_R{read}_trimmed_fastqc.html",
        zip="output/qc/fastqc/{sample}_R{read}_trimmed_fastqc.zip"
    params: ""
    log:
        "output/logs/fastqc/{sample}_R{read}_trimmed.log"
    wrapper:
        "0.34.0/bio/fastqc"

### alignment ###

rule bowtie2:
    input:
        r1="output/trimmed/{sample}_R1.fastq",
        r2="output/trimmed/{sample}_R2.fastq"
    output:
        temp("output/mapped/{sample}.bam")
    conda:
        "envs/bowtie2.yaml"
    params:
        index="/data/bowtie2-indexfiles/genome/S288C_reference_sequence_R64-4-1_20230830",
        extra="-X 2000 --dovetail"
    log:
        "output/logs/bowtie2/{sample}.log"
    threads:
        8
    shell:
        "bowtie2 --threads {threads} -x {params.index} {params.extra} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -Sbh -F4 -f2 -q30 -@ {threads} - > {output}"

rule samtools_sort:
    input:
      	"output/mapped/{sample}.bam"
    output:
        temp("output/mapped/{sample}.sorted.bam")
    log:
        "output/logs/samtools/{sample}.log"
    params:
        "-m 4G"
    threads:
        8
    wrapper:
        "0.35.0/bio/samtools/sort"

### ATAC processing ###

rule remove_duplicates:
    input:
        "output/mapped/{sample}.sorted.bam"
    output:
        bam="output/mapped/{sample}.dedup.bam",
        metrics="output/logs/picard_dedup/{sample}.dedup.metrics.txt"
    log:
        "output/logs/picard_dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "0.35.0/bio/picard/markduplicates"

rule index_dedup:
    input:
        "output/mapped/{sample}.dedup.bam"
    output:
        temp("output/mapped/{sample}.dedup.bai")
    params:
        ""
    wrapper:
        "0.35.0/bio/samtools/index"

rule remove_offset:
    input:
        bam="output/mapped/{sample}.dedup.bam",
        bai="output/mapped/{sample}.dedup.bai"
    output:
        temp("output/mapped/{sample}.dedup.offset.bam")
    log:
        "output/logs/alignmentSieve/{sample}.log"
    conda:
        "envs/deeptools.yaml"
    threads:
        1
    shell:
        "alignmentSieve --numberOfProcessors {threads} --verbose --ATACshift --bam {input.bam} -o {output} 2> {log}"

rule sort_final:
    input:
        "output/mapped/{sample}.dedup.offset.bam"
    output:
        protected("output/mapped/{sample}.final.bam")
    log:
        "output/logs/samtools/{sample}.log"
    params:
        "-m 4G"
    threads:
        8
    wrapper:
        "0.35.0/bio/samtools/sort"

rule index_final:
    input:
        "output/mapped/{sample}.final.bam"
    output:
        protected("output/mapped/{sample}.final.bai")
    params:
        "" # optional params string
    wrapper:
        "0.35.0/bio/samtools/index"

rule samtools_flagstat:
    input:
        "output/mapped/{sample}.final.bam"
    output:
        "output/mapped/{sample}.final.flagstat"
    wrapper:
        "0.35.0/bio/samtools/flagstat"

### subsample ###

rule subsample:
    input:
        "output/mapped/{sample}.final.bam"
    output:
        protected("output/mapped/{sample}.subsample.bam")
    log:
        "output/logs/subsample/{sample}.log"
    conda:
        "envs/sambamba.yaml"
    threads:
        8
    shell:
        """
        nreads=$(samtools view -c {input})
        rate=$(echo "scale=5;10000000/$nreads" | bc)
        sambamba view -f bam -t 5 --subsampling-seed=42 -s $rate {input} | samtools sort -m 4G -@ 8 -T - > {output} 2> {log}
        """

rule index_subsample:
    input:
        "output/mapped/{sample}.subsample.bam"
    output:
        "output/mapped/{sample}.subsample.bai"
    params:
        "" # optional params string
    wrapper:
        "0.35.0/bio/samtools/index"

### BigWig ###

rule bigwig_final:
    input:
        bam="output/mapped/{sample}.final.bam",
        bai="output/mapped/{sample}.final.bai"
    output:
        "output/mapped/{sample}.final.bw"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output}"

rule bigwig_subsample:
    input:
        bam="output/mapped/{sample}.subsample.bam",
        bai="output/mapped/{sample}.subsample.bai"
    output:
        "output/mapped/{sample}.subsample.bw"
    conda:
        "envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output}"

### Peak Calling ###

rule macs2_full:
    input:
        "output/mapped/{sample}.final.bam"
    output:
        xls="output/peaks/{sample}_peaks.xls"
    params:
        sample="{sample}"
    conda:
        "envs/macs2.yaml"
    log:
        "output/logs/macs2/{sample}.log"
    shell:
        "macs2 callpeak -t {input} -f BAM -g 1.0e+9 --outdir output/peaks -n {params.sample} 2> {log}"

rule macs2_subsample:
    input:
        "output/mapped/{sample}.subsample.bam"
    output:
        xls="output/peaks/{sample}_subsample_peaks.xls"
    params:
        sample="{sample}_subsample"
    conda:
        "envs/macs2.yaml"
    log:
        "output/logs/macs2/{sample}.subsample.log"
    shell:
        "macs2 callpeak -t {input} -f BAM -g 1.0e+9 --outdir output/peaks -n {params.sample} 2> {log}"

### Complete QC ###

rule multiqc: 
    input:
        expand("output/qc/fastqc/{sample}_R{read}_trimmed_fastqc.zip",zip,sample=SAMPLES,read=READS),
        expand("output/logs/bowtie2/{sample}.log", sample=SAMPLES),
        expand("output/peaks/{sample}_peaks.xls", sample=SAMPLES),
        expand("output/peaks/{sample}_subsample_peaks.xls", sample=SAMPLES)
    output:
        "output/qc/multiqc.html"
    log:
        "output/logs/multiqc.log"
    params:
        ""
    wrapper:
        "v3.3.3/bio/multiqc"
