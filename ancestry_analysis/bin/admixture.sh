#This script contains the commands to run the clustering analysis with ADMIXTURE Version 1.3.0 for landraces, modern varieties, Zea mays ssp. mexicana and Zea mays ssp. parviglumis

#To reproduce the ancestry analysis the Zea diploperennis samples must be excluded



for K in {1..8}; do
~/bin/admixture --cv Path/2/bedFile/binaryPlinkFileWhithout_Zeadiploperennis.bed $K | tee Path/2/out/directory/out_admixture/Zea_mays.AGPv4.2x_0.8_0.01.NewID_${K}.out;
done

#Mostrar valores de CV
grep -h CV Path/2/out/directory/Zea_mays.AGPv4.2x_0.8_0.01.NewID_*.out > Path/2/out/directory/Zea_mays.AGPv4.2x_0.8_0.01.NewID_K1a8_value_test.txt
