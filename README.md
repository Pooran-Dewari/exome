# exome

### files
location /exports/eddie/scratch/pdewari/exome
```
index36_CCAACA_L002_R1_001.fastq.gz
index36_CCAACA_L002_R2_001.fastq.gz
MD5.txt

```
### download genome and annotation
```
 wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
 wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
```

### index genome
Note: bwa indexing can't use multi-threading; this should take about 45 min.
```
#!/bin/bash
 
#$ -V -cwd
#$ -l h_rt=04:00:00 
#$ -l h_vmem=32G 
#$ -pe sharedmem 1
#$ -P 
 
. /etc/profile.d/modules.sh

 
module load roslin/bwa/0.7.18
module load igmm/apps/samtools/1.20
 
bwa index GRCh38.primary_assembly.genome.fa
samtools faidx GRCh38.primary_assembly.genome.fa
```
### align
```
#!/bin/bash

#$ -V -cwd
#$ -l h_rt=12:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 4
#$ -P 

. /etc/profile.d/modules.sh


module load roslin/bwa/0.7.18
module load igmm/apps/samtools/1.20

bwa mem -t 4 GRCh38.primary_assembly.genome.fa \
/exports/eddie/scratch/pdewari/exome/index36_CCAACA_L002_R1_001.fastq.gz \
/exports/eddie/scratch/pdewari/exome/index36_CCAACA_L002_R2_001.fastq.gz | \
samtools sort -o exome.sorted.bam

```
