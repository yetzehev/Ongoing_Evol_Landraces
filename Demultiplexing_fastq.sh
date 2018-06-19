#!/bin/sh
#This command was run to demultiplex samples from fastq files with GBSX_v1.3 
#30 agosto 2017
for i in C872HANXXX_7_fastq C7U03ANXX_7_fastq CABV0ANXX_1_fastq CABV0ANXX_2_fastq CABV0ANXX_3_fastq
do
java -jar ./GBSX/releases/latest/GBSX_v1.3.jar --Demultiplexer -f1 ../Zmays/$i/data/$i.gz -i ../Zmays/$i/meta/$i.barcodes.txt -gzip true -o ../Zmays/$i/fastq_demultiplexed/
done


  


