#!/bin/bash
#This script contains the command to perform the Variant Discovery
#I  used a modified version of the scrip https://github.com/owensgl/argentina_helianthus/blob/master/make_gcvf_gatk.bash written by Greg Owens

#This command where executed from bin directory
#1) To do a dictionary with picard.jar

java -jar picard.jar CreateSequenceDictionary R= ~/Zmays/Zea_mays.AGPv4.dna.toplevel.fa O=~/Zmays/Zea_mays.AGPv4.dna.toplevel.fa.dict

#To do an index with samtools
./samtools faidx ~/Zmays/Zea_mays.AGPv4.dna.toplevel.fa


#####Then, everything would be ready
#You have to run the next script, remember modify it according your requirements


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


