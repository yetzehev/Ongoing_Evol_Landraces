#This file containst the commands used to perform the gene flow test from MV1 into LR collected after 2000, when the reference group are the Landraces sampled before 1960; and Zea diploperennis is the outgrop

#An example for the bamlist file, the size population file and the populations file is provided in the meta directory

#Model 1.1: LR<1960 / LR>2000 / MV1 /Zp

~/bin/angsd/angsd -doAbbababa2 1 -bam Path/2/bamfiles/list/Model_1.1.bamlist -sizeFile Path/2/sizePopulation/files/Model_1.1_LR1960_2000.size -doCounts 1  -out Path/2/out/directory/Model1.1 -useLast 1 -nThreads 2 -minQ 20 -minMapQ 30

#Model 1.1 1960_2000_MV1_Zp: Command to compute Patterson's D 
Rscript ~/bin/angsd/R/estAvgError.R angsdFile="Path/2/out/directory/Model1.1" out="Path/2/out/directory/Model1.1.Dstats" sizeFile=Path/2/sizePopulation/files/Model_1.1_LR1960_2000.size nameFile=Path/2/Populations/Names/files/Model_1.1_LR1960_2000.size/Model_1.1.popnames
