## Fst in windows
The Fst indexes were calculated for LR an WR, all were compared versus the two modern varieties genetic groups
1) Modern Varieties cluster 1 (MV1) and 2 (MV2)
2) Landraces collected before 1960
3) Landraces collected between 1960-1980
4) Landraces collected after 2000
5) Wild relatives: Zea mays ssp. parviglumis collected between 1978-1984
6) Wild relatives: Zea mays ssp. parviglumis collected after 2000
7) Wild relatives: Zea mays ssp. mexicana collected between 1983-1984
8) Wild relatives: Zea mays ssp. mexicana collected after 2000
The scripts were written by Gregory Owens and are available at https://github.com/owensgl

####Define  populations
This script requires three files to run.
1) A VCF file with the genotypes to evaluate
2) set1.sampsleinfo: A two columns file separated by tab. The first columns are the sample names as these appear in the VCF file,and the second column is the population name.
i.e.
```
"Sample"	"Pop"
PpNa194302_PpNa1943	PpNa1943
PpNa194303_PpNa1943	PpNa1943
PpNa194304_PpNa1943	PpNa1943
PpNa194305_PpNa1943	PpNa1943
PpNa194306_PpNa1943	PpNa1943
ChNa195806_ChNa1958	ChNa1958
ChNa195807_ChNa1958	ChNa1958
ChNa195808_ChNa1958	ChNa1958
ChNa195805_ChNa1958	ChNa1958
ChNa195809_ChNa1958	ChNa1958
TbNa1946a1_TbNa1946a	TbNa1946a
TbNa1946a2_TbNa1946a	TbNa1946a
TbNa1946a3_TbNa1946a	TbNa1946a
TbNa1946a4_TbNa1946a	TbNa1946a
TbNa1946a5_TbNa1946a	TbNa1946a
ZmNa19441b_ZmNa1944	ZmNa1944
ZmNa194402_ZmNa1944	ZmNa1944
BoH1201501_BoH12015	BoH12015
BoH1201502_BoH12015	BoH12015
BoH1201503_BoH12015	BoH12015
BoH1201504_BoH12015	BoH12015
BoH1201505_BoH12015	BoH12015
```

3) set1.popinfo : A two columms file separated by tab file. First column is the population name,and the second column is the number of group.

i.e.

	"Pop"	"Npop"
	PpNa1943	1
	PpNa1943	1
	PpNa1943	1
	PpNa1943	1
	PpNa1943	1
	ChNa1958	1
	ChNa1958	1
	ChNa1958	1
	ChNa1958	1
	ChNa1958	1
	TbNa1946a	1
	TbNa1946a	1
	TbNa1946a	1
	TbNa1946a	1
	TbNa1946a	1
	ZmNa1944	1
	ZmNa1944	1
	BoH12015	2
	BoH12015	2
	BoH12015	2
	BoH12015	2
	BoH12015	2

## Compute Fst per site
To use the scripts an example is:

	cat data.vcf | perl vcf2vertical_bi_basic.pl | perl SNPtable2Fst.pl set1.sampleinfo set1.popinfo > fst.txt


## To average Fst in windows
The window size parameter should be provided in number of pair of bases

	cat out/fst.txt | perl pop_gen/SlidingWindow_onlyfst.pl 1000000 > out/fst_1Mb.txt
