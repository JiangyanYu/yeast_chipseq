docker run -it --rm --name Chip_preprocessing \
 -v /home/yu.j/sciebo/G4/Zuo1ChIPseq2020/output/:/output/ \
 -v /home/yu.j/sciebo/G4/Zuo1ChIPseq2020/fastq/:/data/ \
 -v /home/yu.j/sciebo/G4/Zuo1ChIPseq2020/envs/:/output/envs/ \
 -v /home/yu.j/sciebo/G4/S288C_reference_genome_R64-4-1_20230830/bowtie_index/:/genome \
 jiangyanyu/yeast_chipseq:1.0 /bin/bash
