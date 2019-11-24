# SNP Calling
## 1. Demultiplexing
The samples demultiplexing was performed with GBSx and the SNP Calling was done following the GATK Best Practices algoritm (https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS&p=1)

Note: Check the directory darjeeling_zoology_ubc/bin/bin_SNP_Calling
##1. Fastq demultiplexing with GBS x
The fastq files were demultiplexed with GBSX_v1.3 (https://github.com/GenomicsCoreLeuven/GBSX) using the next command for all  the fastq files.
 - f1: Is fastq file 
 - i: Is the file C7U03ANXX_7_barcodes.txt, it has three columns whitout header:
A) Sample name
B) Barcode
C) Enzyme used to digest the samples
 - gzip: Isthe format of the fastq file
 - o: Is the path for the output file

```
java -jar ../bin/GBSX/releases/latest/GBSX_v1.3.jar --Demultiplexer -f1 ../data/C7U03ANXX_7_fastq.gz -i ../meta/C7U03ANXX_7_barcodes.txt -gzip true -o ../out_lane4/
```
It step should be done in a couple of hours (or less) per lane using  a server as the CONABIO one: using 5 threats .
## 2. Alignment
The  reads aligment was performed using an edited version for the script **align_process.bash**
availaible in the github of Greg Owens: https://github.com/owensgl/argentina_helianthus/blob/master/align_process.bash

This step requieres the installation of:
- Picard (https://sourceforge.net/projects/picard/)
- miniconda
 ```wget -c http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
   ```
- SamTools, I used samtools-1.5
 ```
wget -c http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh. (https://github.com/samtools/samtools)
 ```
- nextgenmap-0.5.3-0 (https://github.com/Cibiv/NextGenMap)
- The **trimmomatic** directory has  a subdirectory called **adapters** that must contain the sequence of the adapters used to do the libraries. The file with the adapter sequences  must be in fasta format, example:
```
>seq1
CWGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
>seq2
ACACTCTTTCCCTACACGACGCTCTTCCGATCTAGGC
>seq3
ACACTCTTTCCCTACACGACGCTCTTCCGATCTGATT
```

Command line to insert the header for each sequence:
```awk '{print ">seq" NR} 1' file```

The next script was used to align the demultiplexed fastq files sequenced in the lane C872HANXX_7, look for this file  in the directory SNP_Calling_v4/bin

```
fastq="/home/rojas/Zmays/C872HANXX_7/fastq"
ngm="/home/rojas/bin/miniconda3/pkgs/nextgenmap-0.5.3-0/bin/ngm"
ref="/home/rojas/Zmays/Zea_mays.AGPv4.dna.toplevel.fa" # It is the reference genome in gz format
bam="/home/rojas/Zmays/C872HANXX_7/bam" #output directory for bam files
picardtools='/home/rojas/bin/picard.jar' #Path to picard executable file
log="/home/rojas/Zmays/C872HANXX_7/log" # output directory for log info
tmp="/home/rojas/Zmays/C872HANXX_7/tmp" # Path for temporal files
samtools='/home/rojas/bin/samtools-1.5/samtools'#Path to samtools
project="Zmays_v4_2017" #Project name
home="/home/rojas/Zmays/C872HANXX_7" # home directory
bin="/home/rojas/bin" #Directory for the executable files
trimmed="/home/rojas/Zmays/C872HANXX_7/trimmed" #Path for temporal file
javarules="-Djava.io.tmpdir=/home/rojas/speedy/tmp" 
ncores="8" #Number opf cores used to run this process
trim="/home/rojas/bin/trimmomatic" 

ls $fastq |grep -v nobar |  grep fastq | grep R1 | sed s/.R1.fastq.gz// | sort -r | uniq > $home/meta/samplelist.txt

while read name
do
if [  ! -f $bam/$name.sort.bam ]; then 
	java -jar $trim/classes/trimmomatic.jar SE -phred33 $fastq/"$name".R1.fastq.gz $trimmed/"$name".R1.trimmed.fastq ILLUMINACLIP:$trim/adapters/adapters.fa:2:30:10:8:T SLIDINGWINDOW:4:15 MINLEN:36
    echo "Aligning paired $name"

    echo "Aligning unpaired $name"
    $ngm -r $ref -q $trimmed/"$name".R1.trimmed.fastq  -o $bam/$name.bam -t $ncores -b --rg-id $name --rg-sm $name --rg-pl illumina --rg-pu $project --rg-lb C872HANXX_7


    echo "Processing for $name"

        $samtools sort -@ $ncores $bam/$name.bam > $bam/${name}.sort.bam 
	$samtools index $bam/${name}.sort.bam
	rm $trimmed/"$name".R1.trimmed.fastq
	rm $bam/$name.bam
}>& $log/$name.alignment.log
fi
done < $home/meta/samplelist.txt
```


## 3. Variant Discovery 
The variant discovery was run with the Haplotype Caller tool of GATK. I used a modified version of the script https://github.com/owensgl/argentina_helianthus/blob/master/make_gcvf_gatk.bash available at Greg Owens github.

This scripts maps the bam fileS to the genome reference using GATK Haplotype Caller
##### WARNINGS:
Before to perform the raw SNP Calling you hava to do 2 previos steps:
#### 1) To create a dictionary with picard.jar
```
java -jar picard.jar CreateSequenceDictionary R= ~/Zmays/Zea_mays.AGPv4.dna.toplevel.fa O=~/Zmays/Zea_mays.AGPv4.dna.toplevel.fa.dict
```
#### 2) To create an index with samtools
```
./samtools faidx ~/Zmays/Zea_mays.AGPv4.dna.toplevel.fa
```
##### Then, everything will be ready
You have to run the next script, remember to modify it according to your requirements

```
#!/bin/sh

javarules="-Djava.io.tmpdir=/home/rojas/speedy/tmp"
ncores="15"

while read prefix
do
 if [ ! -f /home/rojas/Zmays/gvcf/$prefix.gvcf.vcf.idx ]
then

#Call GATK HaplotypeCaller
java -Xmx50g -jar /home/rojas/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
	-R /home/rojas/Zmays/Zea_mays.AGPv4.dna.toplevel.fa \
	-log /home/rojas/Zmays/log/$prefix.HaplotypeCaller.log \
	-T HaplotypeCaller \
	-nct $ncores \
	-I /home/rojas/Zmays/C7U03ANXX_7/bam/$prefix.sort.bam \
	--emitRefConfidence GVCF \
	--max_alternate_alleles 3 \
	-variant_index_type LINEAR \
	-variant_index_parameter 128000 \
	-o /home/rojas/Zmays/gvcf/$prefix.gvcf.vcf ;
fi
done < /home/rojas/Zmays/C7U03ANXX_7/meta/samplelist.txt
```
## 4. Merge gvcf files
The last step performs the SNP Calling for each sample, after this you must merge the gvcf files to get a final VCF file with all the samples and SNPs. 
- Save all the gvcf in the same directory
- I used a modified version of the Greg Owens script available at:
https://github.com/owensgl/argentina_helianthus/blob/master/Genotype_gvcf.sh
- Remember modify it according to your needs

```
combinedGVCFs='/home/rojas/Zmays/gvcf'
ls $combinedGVCFs | grep "vcf" | grep -v ".idx"   > GVCFs.samplelist.txt
tmp=""
while read prefix
do
        tmp="$tmp --variant $combinedGVCFs/$prefix"
done < GVCFs.samplelist.txt

java -Xmx40g -jar /home/rojas/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
        -nt 15 \
        -l INFO \
        -R /home/rojas/Zmays/Zea_mays.AGPv4.dna.toplevel.fa \
        -log /home/rojas/Zmays/log/GenotypeGVCFsAGPv4.log \
        -T GenotypeGVCFs \
        $tmp \
        -o /home/rojas/Zmays/Zea_mays.AGPv4.vcf \
	-hets 0.01 \
        --max_alternate_alleles 4
```

##### The output from this step is ready to be filter!!!
