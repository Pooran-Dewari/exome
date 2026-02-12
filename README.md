### exome-seq analysis


##### files
location /exports/eddie/scratch/pdewari/exome
```
index36_CCAACA_L002_R1_001.fastq.gz
index36_CCAACA_L002_R2_001.fastq.gz
MD5.txt

```


##### download genome and annotation
```
 wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
 wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
```


##### index genome
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



##### align
Takes about 20 min.
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
##### mark duplicates
If picard gives error, try loading the picard module in eddie scratch and check it works, then submit the qsub.  
Takes about 10 min.
```
#!/bin/bash

#$ -V -cwd
#$ -l h_rt=06:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 2
#$ -P 




module load igmm/apps/picard/3.1.1


picard MarkDuplicates \
-I exome.sorted.bam \
-O exome.dedup.bam \
-M metrics.txt \
--CREATE_INDEX true \
--VALIDATION_STRINGENCY SILENT

module load igmm/apps/samtools/1.20

samtools index  exome.dedup.bam

```

##### variant calling

```
#!/bin/bash

#$ -V -cwd
#$ -l h_rt=06:00:00
#$ -l h_vmem=16G
#$ -pe sharedmem 2
#$ -P 
module load roslin/bcftools/1.20
module load roslin/bedtools/2.31.1

# get exons
awk '$3=="exon"' gencode.v49.annotation.gtf | awk '{
  match($0,/gene_name "([^"]+)"/,a);
  if(a[1]=="CPT1A")
    print $1"\t"($4-1)"\t"$5"\t"a[1]
}' > CPT1A.exons.bed

# merge bed files if duplicate entries
bedtools sort -i CPT1A.exons.bed | bedtools merge -i - > CPT1A.exons.merged.bed


# variant calling
bcftools mpileup \
-f GRCh38.primary_assembly.genome.fa \
-R CPT1A.exons.merged.bed \
exome.dedup.bam | \
bcftools call -mv -Oz -o CPT1A.vcf.gz

# index
bcftools index CPT1A.vcf.gz

# generate consensus seq
bcftools consensus \
-f GRCh38.primary_assembly.genome.fa \
-R CPT1A.exons.bed \
CPT1A.vcf.gz > CPT1A.sample.fa
```


##### check variant
```
$ bcftools view CPT1A.vcf.gz | grep -v "^#" 
```
chr11	68781872	.	A	G	225.417	.	DP=98;VDB=0.186076;SGB=-0.693147;MQSBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,41,42;MQ=60	GT:PL	1/1:255,250,0    
Only one single variant, so sampl.fa above could be empty.


##### extract sequence and apply the snp

```
```
