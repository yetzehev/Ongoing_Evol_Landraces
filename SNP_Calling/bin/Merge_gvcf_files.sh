# After run Haplotype Caller you must merge the gvcf files to get a final VCF file with all the samples and SNPs. 
# Save all the gvcf in the same directory
# I used a modified version of the Greg Owens script available at:
# https://github.com/owensgl/argentina_helianthus/blob/master/Genotype_gvcf.sh


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

#The output is ready to be filter
