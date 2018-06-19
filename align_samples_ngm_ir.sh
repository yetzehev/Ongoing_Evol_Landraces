#This script is a modified version of the script available at https://github.com/owensgl/argentina_helianthus/blob/master/align_process.bash
#written by Greg Owens 
fastq="/home/rojas/Zmays/C7U03ANXX_7/fastq"
ngm="/home/rojas/bin/miniconda3/pkgs/nextgenmap-0.5.3-0/bin/ngm"
ref="/home/rojas/Zmays/Zea_mays.AGPv4.dna.toplevel.fa.gz" # Reference genome can be a gz file
bam="/home/rojas/Zmays/C7U03ANXX_7/bam" #bam files pathway
picardtools='/home/rojas/bin/picard.jar'
log="/home/rojas/Zmays/C7U03ANXX_7/log" # output directory
tmp="/home/rojas/Zmays/C7U03ANXX_7/tmp" # output directory
samtools='/home/rojas/bin/samtools-1.5/samtools'
project="Zmays_v4_2017"
home="/home/rojas/Zmays/C7U03ANXX_7" # output directory
bin="/home/rojas/bin"
trimmed="/home/rojas/Zmays/C7U03ANXX_7/trimmed"
javarules="-Djava.io.tmpdir=/home/rojas/speedy/tmp" 
ncores="8"
trim="/home/rojas/bin/trimmomatic"
#This script was modified because I have single strand sequencing data
ls $fastq |grep -v nobar |  grep fastq | grep R1 | sed s/.R1.fastq.gz// | sort -r | uniq > $home/meta/samplelist.txt

while read name
do
if [  ! -f $bam/$name.sort.bam ]; then #adaptors must be in fasta format
{
	java -jar $trim/classes/trimmomatic.jar SE -phred33 $fastq/"$name".R1.fastq.gz $trimmed/"$name".R1.trimmed.fastq ILLUMINACLIP:$trim/adapters/adapters.fa:2:30:10:8:T SLIDINGWINDOW:4:15 MINLEN:36
    echo "Aligning paired $name"

    echo "Aligning unpaired $name"
    $ngm -r $ref -q $trimmed/"$name".R1.trimmed.fastq  -o $bam/$name.bam -t $ncores -b --rg-id $name --rg-sm $name --rg-pl illumina --rg-pu $project --rg-lb C7U03ANXX_7


    echo "Processing for $name"

        $samtools sort -@ $ncores $bam/$name.bam > $bam/${name}.sort.bam 
	$samtools index $bam/${name}.sort.bam
	rm $trimmed/"$name".R1.trimmed.fastq
	rm $bam/$name.bam
}>& $log/$name.alignment.log
fi
done < $home/meta/samplelist.txt

