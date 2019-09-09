######################################################################################
# PCA_Romay_Data (Romay_Rojas_):
#Note: I merged the MV collected in simpatry with LRs and the public data from Romay et al., 2013
#The genome version is v2
########################################################################################
##### Para usar con SNPRelate
#Load libraries
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(permute)
library(lattice)
library(vegan)
library(RColorBrewer)
#To make a gds file from  binary plink files
bed.fn <- ("~/hubic/MisScripts/21.Tree_NJ/data/PCA_P_tree_NID.bed")
fam.fn <- ("~/hubic/MisScripts/21.Tree_NJ/data/PCA_P_tree_NID.fam")
bim.fn <- ("~/hubic/MisScripts/21.Tree_NJ/data/PCA_P_tree_NID.bim")
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn,"~/hubic/MisScripts/21.Tree_NJ/data/PCA_P_tree_NID.gds")
#summary
snpgdsSummary("~/hubic/MisScripts/21.Tree_NJ/data/PCA_P_tree_NID.gds")
# To load the gds file
data <- snpgdsOpen("~/hubic/MisScripts/21.Tree_NJ/data/PCA_P_tree_NID.gds", allow.duplicate = TRUE)
# Check snp.ids
head(read.gdsn(index.gdsn(data, "snp.id")))
# Check sample.ids
head(read.gdsn(index.gdsn(data, "sample.id")))
# To get IDs list
sample.id <- read.gdsn(index.gdsn(data, "sample.id"))
#Make a data frame and rename the column
sample.id <- as.data.frame(sample.id)
colnames(sample.id) <- "Samples"
length(sample.id)
#Populations 
Groups <- read.delim("~/hubic/MisScripts/21.Tree_NJ/meta/PCA_P_breedingPrograms.txt", header = T)
dim(Groups)
GeneticPool.2 = levels(as.factor(Groups$GeneticPool.2))
GeneticPool.1 = levels(as.factor(Groups$GeneticPool.1))
Pop.structure = levels(as.factor(Groups$Pop.structure))
Groups$NSample <- (1:nrow(Groups)) #change this for every data set
dim(Groups)

# PCA
pca <- snpgdsPCA(data, num.thread=2,verbose=TRUE)

# Compute % of variance for EGV-1 and EGV-2
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

x<-round(pc.percent, 2)
sum(x[1:4])
sum(x[1:10])
sum(x[1:30])
#Make a data frame
tab1 <- data.frame(sample.id = Groups$New.ID,
                  Groups=unlist(Groups$GeneticPool.1),
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

#Define pallete
levels(Groups$GeneticPool.1)
levels(Groups$GeneticPool.2)
levels(Groups$Pop.structure)

#Define palletes
b73 <-brewer.pal(5,"Greys")
pools <-brewer.pal(5,"Set2")
mvs.lr <- c("dodgerblue","red","darkgoldenrod1")
Niger.parvi <- brewer.pal(4,"Greens")
pools2 <- brewer.pal(7,"PuOr") 


mycols <- c(b73[2:5], mvs.lr, pools, Niger.parvi[2:4],pools2)
 
#Plot the two first eigenvectors
#png(filename = "~/hubic/MisScripts/7.PCA/out/plots/PCA_P_geneticPool_1.png", height = 500, width = 600)
#Color by group
p1 <- ggplot(tab1, aes(x=EV1, y=EV2)) +
  geom_point(aes(colour=Groups),size =1.5) + 
  scale_color_manual(values=mycols,
                     breaks = c("Landraces","MV1","MV2","Tropical",
                                "B37","B73","B84","Costeño",    
                                "Non-stiff stalk","parviglumis","popcorn","stiff stalk",
                                "sweet corn", "Tuxpeño","Zapalote chico" ))+
  xlab(paste0("eigenvector 1 (", round(pc.percent, 2)[1], "% )")) +
  ylab(paste0("eigenvector 2 (", round(pc.percent, 2)[2], "% )"))+
  theme(axis.text.x = element_text(colour="grey20",size=12),
        axis.text.y = element_text(colour="grey20",size=12),  
        axis.title.x = element_text(colour="grey20",size=15),
        axis.title.y = element_text(colour="grey20",size=15),
        legend.title =element_text(size=12),
        legend.text=element_text(size=12)) +
  theme(plot.title = element_text(size = 10))+
  labs (title = "PCA color by Genetic Pool (1002 samples, 13953 SNPs)")
p1
#dev.off()
