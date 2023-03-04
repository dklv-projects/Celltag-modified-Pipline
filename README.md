# Celltag-modified-pipline使用说明
测试文件下载地址为https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR7347033&display=metadata
#下载文件，文件大小为16.6G
wget -c https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR7347033/SRR7347033
#下载文件为bam格式，需转为fastq
#通过conda安装bamfastq
conda install -c bioconda 10x_bamtofastq
#安装后pwd查看目录，运行bamtofastq
~/anaconda3/envs/xx/bin/bamtofastq "/home/visitor/anaconda3/bin/hf1.d15.possorted_genome_bam.bam.1" /home/visitor/anaconda3/bin/fastq/
#cellranger count命令前，需提前构建比对基因组
#下载并解压mouse gtf
wget -c https://ftp.ensembl.org/pub/release-108/gtf/mus_musculus/Mus_musculus.GRCm39.108.gtf.gz
gunzip Mus_musculus.GRCm39.108.gtf.gz
#下载并解压mouse fastq
wget -c https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
#新建fa文件并添加celltag
vim CellTag.UTR #vim进去后，I编辑 esc退出编辑 ：wq保存后退出）
>CellTag.UTR
GAATTCGATGACAGGCGCAGCTTCCGAGGGATTTGAGATCCAGACATGATAAGATACATT
GATGAGTTTGGACAAACCAAAACTAGAATGCAGTGAAAAAAATGCCTTATTTGTGAAATT
TGTGATGCTATTGCCTTATTTGTAACCATTATAAGCTGCAATAAACAAGTTAACA
>GFP.CDS
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCTGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGCGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA
cat Mus_musculus.GRCm39.dna.primary_assembly.fa CellTag.fa >mouse.celltag.fa
#新建gtf文件并添加
vim GTF 
#GTF Entries:
CellTag.UTR     custom     exon     1     175     .     +     .     gene_id "CellTag.UTR"; transcript_id "celltag.utr";
GFP.CDS     custom     exon     1     720     .     +     .     gene_id "GFP.CDS"; transcript_id "gfp.cds”;）
cat Mus_musculus.GRCm39.108.gtf GTF >mouse.celltag.gtf
#构建比对基因组
cellranger mkref \--genome=Mus_musculus_genome \--fasta=mouse.celltag.fa \--genes=mouse.celltag.gtf
#cellranger count
cellranger count  --id=sampleT1 \
--transcriptome=/home/visitor/anaconda3/bin/Mus_musculus_genome  \
--fastqs=/home/visitor/anaconda3/bin/fastq \
--sample=bamtofastq \
--localcores=8 \
--localmem=64 \
#选择cellranger outs中的possorted_genome_bam.bam和barcode.tsv文件进行下一步分析
#Bam文件是二进制文件，需要先安装samtools（先安装brew，再用brew安装samtools）
/bin/zsh -c "$(curl -fsSL https://gitee.com/cunkai/HomebrewCN/raw/master/Homebrew.sh)” 
brew install samtools
#安装grep,grep运行前需配置环境
PATH="/opt/homebrew/opt/grep/libexec/gnubin:$PATH"
#bash
#MEF
samtools view possorted_genome_bam.bam | grep -P 'GGT[ACTG]{8}GAATTC' > v1.celltag.reads.out
#D3
samtools view possorted_genome_bam.bam | grep -P 'GTGATG[ACTG]{8}GAATTC' > v2.celltag.reads.out
#D13
samtools view possorted_genome_bam.bam | grep -P 'TGTACG[ACTG]{8}GAATTC' > v3.celltag.reads.out
#选用上一步的output celltag.reads.out进行下一步分析
#gawk也需要配置环境
 PATH="/opt/homebrew/opt/gawk/libexec/gnubin:$PATH"
#bash
#MEF
./scripts/celltag.parse.reads.10x.sh -v tagregex="CCGGT([ACTG]{8})GAATTC" v1.celltag.reads.out > v1.celltag.parsed.tsv
#D3
./scripts/celltag.parse.reads.10x.sh -v tagregex="GTGATG([ACTG]{8})GAATTC" v2.celltag.reads.out > v2.celltag.parsed.tsv
#D13
./scripts/celltag.parse.reads.10x.sh -v tagregex="TGTACG([ACTG]{8})GAATTC" v3.celltag.reads.out > v3.celltag.parsed.tsv
#选用上一步的output celltag.parsed.tsv和cellranger outs的barcode.tsv进行下一步分析
#bash
#MEF
Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v1.celltag.parsed.tsv hf1.d15.v1
#D3
Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v2.celltag.parsed.tsv hf1.d15.v2
#D13
Rscript ./scripts/matrix.count.celltags.R ./cell.barcodes/hf1.d15.barcodes.tsv v3.celltag.parsed.tsv hf1.d15.v3
#剩余的分析在R中进行
#clonecalling
#安装需要的R包
install.packages("igraph")
install.packages("proxy")
install.packages("corrplot")
install.packages("data.table")
#调用所需要的包
library(igraph)
library(proxy)
library(corrplot)
library(data.table)
source("./scripts/CellTagCloneCalling_Function.R")
#运行clonecalling.R脚本
#Lineage and Network Visualization
#安装需要的R包
install.packages("tidyverse")
install.packages("foreach")
install.packages("networkD3")
#调用需要的R包
library(tidyverse)
library(foreach)
library(networkD3)
source("./scripts/function_source_for_network_visualization.R")
source("./scripts/function_source_for_network_construction.R")
#运行linagetree.R脚本
