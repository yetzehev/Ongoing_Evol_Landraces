######################################################################################
# PCA_Romay_Data (Romay_Rojas_):
#Note: I merged the IL that I collected and the public data from CIMMYT
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
#To make a gds file from a binary plink file
bed.fn <- ("~/hubic/MisScripts/7.PCA/data/PCA_0_Rojas_Romay_CMLs_LRs.bed")
fam.fn <- ("~/hubic/MisScripts/7.PCA/data/PCA_0_Rojas_Romay_CMLs_LRs.fam")
bim.fn <- ("~/hubic/MisScripts/7.PCA/data/PCA_0_Rojas_Romay_CMLs_LRs.bim")
# convert
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn,"~/hubic/MisScripts/7.PCA/data/PCA_0_Rojas_Romay_CMLs_LRs.gds")
#summary
snpgdsSummary("~/hubic/MisScripts/7.PCA/data/PCA_0_Rojas_Romay_CMLs_LRs.gds")
# To load gds file
data <- snpgdsOpen("~/hubic/MisScripts/7.PCA/data/PCA_0_Rojas_Romay_CMLs_LRs.gds", allow.duplicate = TRUE)
# Check snp.ids
head(read.gdsn(index.gdsn(data, "snp.id")))
# Check sample.ids
head(read.gdsn(index.gdsn(data, "sample.id")))
#To get IDs list
sample.id <- read.gdsn(index.gdsn(data, "sample.id"))
#Make a data frame and rename the column
sample.id <- as.data.frame(sample.id)
colnames(sample.id) <- "Samples"
length(sample.id)
#Populations 
Groups <- read.delim("~/hubic/MisScripts/21.Tree_NJ/meta/PCA.O.maize.classification.txt",header = T)
dim(Groups)
levels(as.factor(Groups$GeneticPool.2))
Groups$NSample <- (1:nrow(Groups)) #change this for every data set
dim(Groups)

# PCA
pca <- snpgdsPCA(data, num.thread=2,verbose=TRUE)

# To compute the % of variance for EGV-1 and EGV-2
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

x<-round(pc.percent, 2)
sum(x[1:4])
sum(x[1:10])
sum(x[1:30])
#Make a data frame
tab <- data.frame(sample.id = Groups$New.ID,
                  Groups=unlist(Groups$GeneticPool.1),
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

#Define pallete
levels(Groups$GeneticPool.2)
levels(Groups$GeneticPool.1)
levels(Groups$Classification.1)
levels(Groups$Classification.3)

 
#Set Color using as reference GeneticPool.1
 b73 <- brewer.pal(7,"Greys")
 lr.purple <- brewer.pal(4,"PRGn") 
 mvs.lr <- c("dodgerblue","red","darkgoldenrod1")
 pools <-brewer.pal(7,"Accent")

  mycols <- c(b73[5:7], lr.purple[1], mvs.lr, pools[1:3],pools[7],pools[5], 
              lr.purple[2], "bisque3",lr.purple[3:4])

#Plot the two first eigenvectors
#png(filename = "~/hubic/MisScripts/7.PCA/out/plots/PCA_O_geneticPool_1.png", height = 500, width = 600)
#Color by group
p <- ggplot(tab, aes(x=EV1, y=EV2)) +
  geom_point(aes(colour=Groups),size =2) + 
  scale_color_manual(values=mycols,
                     breaks =c(
                               "MV1","MV2","Landraces","Tropical",
                               "Costeño", "Tuxpeño","Zapalote chico",
                               "B37","B73","B84", 
                               "Non-stiff stalk","Popcorn","Stiff stalk",
                               "Sweet corn", "Unclassified"
                               ))+
  xlab(paste0("eigenvector 1 (", round(pc.percent, 2)[1], "% )")) +
  ylab(paste0("eigenvector 2 (", round(pc.percent, 2)[2], "% )"))+
  theme(axis.text.x = element_text(colour="grey20",size=12),
        axis.text.y = element_text(colour="grey20",size=12),  
        axis.title.x = element_text(colour="grey20",size=15),
        axis.title.y = element_text(colour="grey20",size=15),
        legend.title =element_text(size=12),
        legend.text=element_text(size=12))+
  ggtitle("PCA_O 2,579 samples, 13,953 SNPs" )
p
#dev.off()

