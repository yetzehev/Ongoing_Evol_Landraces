#!/bin/bash

######################################
#	Step 1: Finding a 'global estimate' of the SFS
#        -doSaf          0
# 	1: perform multisample GL estimation
#	GLs genotype likelihoods
#	1) The recalibrating SOAPsnp model 
#	2) The original GATK model
#	3) SAMtools 1.16+ modified Maq model 
#	4) The type specific error model


#analysisEstLikes.cpp:
	#-GL=0: 
	#1: SAMtools ****
	#2: GATK
	#3: SOAPsnp
	#4: SYK
# I used -GL 1 Samtools model

~/bin/angsd/angsd -bam /home/idalia/hubic/MisScripts/4.Diversity_Pi/meta/LR1960.angsd.bamlist -doSaf 1 -GL 1 -anc ~/hubic/MisScripts/4.Diversity_Pi/meta/Zea_diploperennis.fa.gz -out /home/idalia/hubic/MisScripts/4.Diversity_Pi/out/LR1960.SFS

#Make a loop for all the samples
for  i in MV1 MV2 LR1980 Zmx1980 Zpr1980 Zmx2000 Zpr2000 LR2000; do
~/bin/angsd/angsd -bam /home/idalia/hubic/MisScripts/4.Diversity_Pi/meta/$i.angsd.bamlist -doSaf 1 -GL 1 -anc ~/hubic/MisScripts/4.Diversity_Pi/meta/Zea_diploperennis.fa.gz -out /home/idalia/hubic/MisScripts/4.Diversity_Pi/out/$i.SFS;
done


#- The output from this first step are the file file.saf.pos.gz; file.saf.idx; file.saf.gz
#I ran this step in ayocote9, it used 15g of RAM
for  i in MV1 MV2 LR1980 Zmx1980 Zpr1980 Zmx2000 Zpr2000 LR2000; do
~/bin/angsd/misc/realSFS $i.SFS.saf.idx -P 4 > $i.sfs;
done


######################################
#	Step 2: Calculate the thetas for each site

#	-pest (null) (prior SFS) output from step 1


~/bin/angsd/angsd -bam /home/idalia/hubic/MisScripts/4.Diversity_Pi/meta/LR1960.angsd.bamlist  -doThetas 1 -doSaf 1 -pest /home/idalia/hubic/MisScripts/4.Diversity_Pi/out/LR1960.SFS -anc ~/hubic/MisScripts/4.Diversity_Pi/meta/Zea_diploperennis.fa.gz -GL 1 -out LR1960.theta_site


~/bin/angsd/misc/thetaStat print ~/hubic/MisScripts/4.Diversity_Pi/out/LR1960.theta_site.thetas.idx 2>/dev/null |head #Muestra el encabezado del archivo

######################################
#Step 3a: Estimate Tajimas D and other statistics (window=10Kb; step 1kb)

~/bin/angsd/misc/thetaStat do_stat ~/hubic/MisScripts/4.Diversity_Pi/out/LR1960.theta_site.thetas.idx

~/bin/angsd/misc/thetaStat do_stat ~/hubic/MisScripts/4.Diversity_Pi/out/LR1960.theta_site.thetas.idx -win 5000 -step 1000  -outnames LR1960.thetasWindow.gz

#Step 3b: Estimate Tajimas D and other statistics, I used bigger windows because  I coul not plot the last set (window=100Kb; step 25kb)


~/bin/angsd/misc/thetaStat do_stat /home/idalia/Nuc_Diversity/out/theta_site/LR1960.theta_site.thetas.idx -win 100000 -step 25000  -outnames ./window/LR1960.100kb_25kb

# Step 3c, try another windows to test if the result is associated to the window and step size (window=50Kb; step 10kb)4
#R can habdle this window size

~/bin/angsd/misc/thetaStat do_stat /home/idalia/Nuc_Diversity/out/theta_site/LR1960.theta_site.thetas.idx -win 50000 -step 10000  -outnames /home/idalia/Nuc_Diversity/out/theta_site/window/LR1960.50kb_10kb

# Step 3d, try another windows to test if the result is associated to the window and step size (window=20kb; step 5kb)4
#R can habdle this window size
for LR in LR1960 LR1980 LR2000; do
~/bin/angsd/misc/thetaStat do_stat /home/idalia/Nuc_Diversity/out/theta_site/$LR.theta_site.thetas.idx -win 20000 -step 5000  -outnames /home/idalia/Nuc_Diversity/out/theta_site/window/$LR.20kb_5kb;
done

