######################################################################################
# PCA_Romay_Data (Romay_Rojas_):
#Note: I merged the MV collected in simpatry with LRs and the public data from Romay et al., 2013
#The genome version is v2
########################################################################################
#Load libraries
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(permute)
library(lattice)
library(vegan)
library(RColorBrewer)
#To make a gds file froom binary plink files
bed.fn <- ("~/hubic/MisScripts/7.PCA/data/PCA_Q_Tropical.bed")
fam.fn <- ("~/hubic/MisScripts/7.PCA/data/PCA_Q_Tropical.fam")
bim.fn <- ("~/hubic/MisScripts/7.PCA/data/PCA_Q_Tropical.bim")
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn,"~/hubic/MisScripts/7.PCA/data/PCA_Q_Tropical.gds")
#summary
snpgdsSummary("~/hubic/MisScripts/7.PCA/data/PCA_Q_Tropical.gds")
# To load file
data <- snpgdsOpen("~/hubic/MisScripts/7.PCA/data/PCA_Q_Tropical.gds", allow.duplicate = TRUE)
# Check snp.ids
head(read.gdsn(index.gdsn(data, "snp.id")))
# Check sample.ids
head(read.gdsn(index.gdsn(data, "sample.id")))
# To make an ID list
sample.id <- read.gdsn(index.gdsn(data, "sample.id"))
#Make a data frame and rename the column
sample.id <- as.data.frame(sample.id)
colnames(sample.id) <- "Samples"
length(sample.id)
#Populations 
Groups <- read.delim("~/hubic/MisScripts/7.PCA/meta/PCA_Q_tropicalBreedingPrograms.txt", header = T)
dim(Groups)
GeneticPool.2 = levels(as.factor(Groups$GeneticPool.2))
GeneticPool.1 = levels(as.factor(Groups$GeneticPool.1))
Breeding.program = levels(as.factor(Groups$Breeding.program))
Groups$NSample <- (1:nrow(Groups)) #change this for every data set
dim(Groups)

# PCA
pca <- snpgdsPCA(data, num.thread=2,verbose=TRUE)

# To compute the % of variance for the EGV 1 and EGV 2
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

x<-round(pc.percent, 2)
sum(x[1:4])
sum(x[1:10])
sum(x[1:30])
#Make a data frame
#Omit parviglumis
Groups <-subset(Groups, Breeding.program != "parviglumis")
dim(Groups)
tab1 <- data.frame(sample.id = Groups$New.ID,
                  Groups=unlist(Groups$Breeding.program),
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

#Define pallete
levels(Groups$GeneticPool.1)
levels(Groups$GeneticPool.2)
levels(Groups$Breeding.program)

#Define palletes
mycols <- c("chocolate4","dodgerblue","red","darkgoldenrod1","deeppink3","gray")

#Plot the two first eigenvectors
#png(filename = "~/hubic/MisScripts/7.PCA/out/plots/PCA_Q_breedingProgram.png", height = 500, width = 600)
#Color by group
p1 <- ggplot(tab1, aes(x=EV1, y=EV2)) +
  geom_point(aes(colour=Groups),size =2) + 
  scale_color_manual(values=mycols,
                     breaks = c("MV1","MV2",
                                "Mexico","Nigeria","Cameroon",
                                "Other"))+
  xlab(paste0("eigenvector 1 (", round(pc.percent, 2)[1], "% )")) +
  ylab(paste0("eigenvector 2 (", round(pc.percent, 2)[2], "% )"))+
  theme(axis.text.x = element_text(colour="grey20",size=12),
        axis.text.y = element_text(colour="grey20",size=12),  
        axis.title.x = element_text(colour="grey20",size=15),
        axis.title.y = element_text(colour="grey20",size=15),
        legend.title =element_text(size=12),
        legend.text=element_text(size=12)) +
  labs (title = "PCA_Q Tropical lines, color by Genetic Pool (463 samples, 13953 SNPs)")
p1
#dev.off()

