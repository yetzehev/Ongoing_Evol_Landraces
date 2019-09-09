#2.06.2019
#This script was used to run an admixture analysis with ADMIXTURE version 1.3.0 for 1002 genotypes which included the genotypes from the tropical breeding programs and the MVs collected in simpatry with LRs

for K in {1..20}; do
~/bin/admixture --cv Path/2/bedFile/data/PCA.Tropical.MVs.bed $K | tee Path/2/out/directory/MVs_tropicalLines_Romay_Rojas_${K}.out;
done

#Mostrar valores de CV
grep -h CV Path/2/out/directory/MVs_tropicalLines_Romay_Rojas_*.out > Path/2/out/directory/MVs_tropicalLines_Romay_Rojas_K1a20_value_test.txt
