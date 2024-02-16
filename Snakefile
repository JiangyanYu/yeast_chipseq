"""
Author: Jiangyan Yu, modified from jsschrepping/bioinfo-base-image:jss_v0.0.3
Affiliation: LIMES, Uni Bonn
Aim: A simple Snakemake workflow to process single-end stranded Chip-Seq.
Date: 20240103
Protocol: Followed a protocol from hbctraining (https://github.com/hbctraining/Intro-to-ChIPseq/blob/master/lessons/05_peak_calling_macs.md)
Docker: docker run -it --rm --name Chip_preprocessing \
 -v /home/yu.j/sciebo/G4/Zuo1ChIPseq2020/output/:/output/ \
 -v /home/yu.j/sciebo/G4/Zuo1ChIPseq2020/fastq/:/data/ \
 -v /home/yu.j/sciebo/G4/Zuo1ChIPseq2020/envs/:/output/envs/ \
 -v /home/yu.j/sciebo/G4/S288C_reference_genome_R64-4-1_20230830/bowtie_index/:/genome \
 jiangyanyu/yeast_chipseq:1.0 /bin/bash
Torun: snakemake --snakefile /output/Snakefile --use-conda --cores 28 --verbose
# bowtie2 index basename is the directory plus file name before .1.bt2
"""

#SAMPLES=["SRR11625406_IP_replicate","SRR11625407_IP_replicate","SRR11625408_input_replicate","SRR11625409_input_replicate"]
#IP = ["SRR11625406","SRR11625407"]
#INPUT = ["SRR11625408","SRR11625409"]

SAMPLES = ["SRR11625406_converted_IP_replicate","SRR11625407_converted_IP_replicate","SRR11625408_converted_input_replicate","SRR11625409_converted_input_replicate"]
IP = ["SRR11625406_converted_IP_replicate","SRR11625407_converted_IP_replicate"]
INPUT = ["SRR11625408_converted_input_replicate","SRR11625409_converted_input_replicate"]

READS=["1"] # it is crucial to recognize sample and reads information for snakemake wildcard#

### Target rule ###

rule all:
    input:
        expand("output/qc/fastqc/{sample}_R{read}_fastqc.html", sample=SAMPLES, read=READS),
        expand("output/mapped/{sample}.bam", sample=SAMPLES),
        expand("output/mapped/{sample}.sorted.bam", sample=SAMPLES),
        expand("output/mapped/{sample}.filtered.bam", sample=SAMPLES),
        expand("output/mapped/{sample}.filtered.bai", sample=SAMPLES),
        expand("output/peaks/IP{ip}_INPUT{input}_peaks.xls",ip = IP, input = INPUT)
        #"output/qc/multiqc.html"

### raw fastq ###

rule fastqc_raw:
    input:
        "/data/{sample}_R{read}.fastq"
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


### alignment ###

rule bowtie2:
    input:
        r1="/data/{sample}_R1.fastq"
        #r2="output/trimmed/{sample}_R2.fastq"
    output:
        "output/mapped/{sample}.bam"
    conda:
        "envs/bowtie2.yaml"
    params:
        index="/genome/S288C_reference_sequence_R64-4-1_20230830",
        extra="-X 2000 --dovetail"
    log:
        "output/logs/bowtie2/{sample}.log"
    threads:
        26
    shell:
        #"bowtie2 --threads {threads} -x {params.index} {params.extra} -1 {input.r1} -2 {input.r2} 2> {log} | samtools view -Sbh -F4 -f2 -q30 -@ {threads} - > {output}"
        #"bowtie2 --threads {threads} -x {params.index} {params.extra} -U {input.r1}  2> {log} | samtools view -Sbh -F4 -f2 -q30 -@ {threads} - > {output}"
        "bowtie2 --threads {threads} -x {params.index} {params.extra} -U {input.r1}  2> {log} | samtools view -Sbh -@ {threads} - > {output}"

rule samtools_sort:
    input:
      	"output/mapped/{sample}.bam"
    output:
        "output/mapped/{sample}.sorted.bam"
    log:
        "output/logs/samtools/{sample}.log"
    params:
        "-m 4G"
    threads:
        8
    wrapper:
        "0.35.0/bio/samtools/sort"
        
rule sambamba_filter:
    input:
        "output/mapped/{sample}.sorted.bam"
    output:
        "output/mapped/{sample}.filtered.bam"
    log:
        "output/logs/mapped/{sample}.filter.log"
    conda:
        "envs/sambamba.yaml"
    threads:
        8
    shell:
        """
        sambamba view -h -t 6 -f bam -F "[XS] == null and not unmapped  and not duplicate " {input}  > {output} 2> {log}
        """
        
rule index_align:
    input:
        "output/mapped/{sample}.filtered.bam"
    output:
        "output/mapped/{sample}.filtered.bai"
    params:
        ""
    wrapper:
        "0.35.0/bio/samtools/index" 
        
### Peak Calling ###

rule macs2_full:
    input:
        ip = "output/mapped/{ip}.filtered.bam",
        input = "output/mapped/{input}.filtered.bam"
    output:
        xls="output/peaks/IP{ip}_INPUT{input}_peaks.xls"
    params:
        sample="IP{ip}_INPUT{input}"
    conda:
        "envs/macs2.yaml"
    log:
        "output/logs/macs2/IP{ip}_INPUT{input}.log"
    shell:
        "macs2 callpeak -t {input.ip} -c {input.input} -f BAM -g 1.0e+9 --outdir output/peaks -n {params.sample} "           

