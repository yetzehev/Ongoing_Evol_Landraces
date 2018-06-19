#Load libraries
library(psych)
library(purrr)
library(data.table)
library(gdata)
library(ggplot2)
library(dplyr)
library(tidyr)
library(easyGgplot2)
setwd("~/hubic/MisScripts/15.Fst_Windows/out/RelabeledSamples_15022018/windows_5kb/")
#Load 
file.list <- list.files(pattern='*.5kb.out') #Make a list of files with the pattern that yopu choose
file.list <- setNames(file.list, file.list) # only needed when you need an id-column with the file-names; from purrr package
#Read all the files at once          
Fst_5kb <- map_df(file.list, read.delim, .id = "id") #From purr package, it is similar to lapply
                                                    #Transform their input by applying a function to each element
Fst_5kb$Fst[is.na(Fst_5kb$Fst)] <- 0 # Substitute NA with zeros

###################
by(Fst_5kb[, 8], Fst_5kb[,"id"], summary) #Make summary by pop id
SamplingUnit.Summary <- capture.output(by(Fst_5kb[, 8], Fst_5kb[,"id"], describe))
write.table(SamplingUnit.Summary, "~/hubic/MisScripts/15.Fst_Windows/out/Stats/SamplingUnit.RelabeledSamples.Summary",sep = "\t")
Fst_5kb$id <-as.character(unlist(strsplit(Fst_5kb$id, '\\.5kb.out'))) #Remove the file extension
colnames(Fst_5kb)


################################################################################################
#         STATISTICS
################################################################################################
### #### Landraces: Chalqueño
Fst_ChNa <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_1_ChNa", ] #Select the columns with id and Fst value
class(Fst_ChNa$id)
lr.pops <- levels(as.factor(Fst_ChNa$id)) #Popnames
lr.pops
#Perform a Mann-Whitney test for Chalqueño landrace collected at different periods
Fst_1_ChNa.stat <- kruskal.test(x=Fst_ChNa$Fst,g= as.factor(Fst_ChNa[,1]))
Fst_1_ChNa.RelabeledSamples.ph <-capture.output(pairwise.wilcox.test(Fst_ChNa$Fst,
                                  Fst_ChNa$id,
                                  p.adjust.method="bonferroni"),file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_1_ChNa.RelabeledSamples.ph")
####  Boxplot Chalqueño
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/ChNaPeriod_bp.png", width = 700, height = 300)
fill.lr <- c ("dodgerblue4","dodgerblue4","dodgerblue","dodgerblue","blueviolet","blueviolet","firebrick3","firebrick3")
line.lr <- c ("gray57")
Fst_ChNa_bp <- ggplot(Fst_ChNa, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)

Fst_ChNa_bp + scale_x_discrete(labels=c("ChNa1958-IL1","ChNa1958-IL2", 
                                        "ChNa1972-IL1", "ChNa1972-IL2",
                                        "ChNa2003-IL1","ChNa2003-IL2",
                                        "ChNa2015-IL1","ChNa2015-IL2"))+
 scale_fill_manual(values=fill.lr)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
        size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/ChNaPeriod_vp.png", width = 700, height = 300)
a <- ggplot2.violinplot(data=Fst_ChNa, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="black", meanPointFill="black",
                        groupName = 'id',
                        groupColors = fill.lr)+guides(fill=FALSE)+
theme(text = element_text(size=6),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=11),
      axis.line = element_line(colour = "black", 
                               size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(labels=c("ChNa1958-IL1","ChNa1958-IL2", "ChNa1972-IL1", "ChNa1972-IL2",
                            "ChNa2003-IL1","ChNa2003-IL2","ChNa2015-IL1","ChNa2015-IL2"))
a
dev.off()

####  #### WR in sympatry with Chalqueño
Fst_ChTm <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_1_ChTm", ]
wr.pops <- levels(as.factor(Fst_ChTm$id)) #Popnames
Fst_1_ChTm.stat <- kruskal.test(x=Fst_ChTm$Fst,g= as.factor(Fst_ChTm[,1]))
capture.output(pairwise.wilcox.test(Fst_ChTm$Fst,Fst_ChTm$id,
                                                    p.adjust.method="bonferroni"),
               file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_1_ChTm.RelabeledSamples.ph")
####Boxplot WR in sympatry with Chalqueño
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/ChTmPeriod.png", width = 600, height = 300)
fill_mx <- c("red","red","black ","black","dimgray","dimgray")
line_mx <- "gray57"
Fst_ChTm_bp <- ggplot(Fst_ChTm, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)
Fst_ChTm_bp + scale_x_discrete(labels=c("ChTm1984-IL1","ChTm1984-IL2","ChTm2003-IL1",
                                        "ChTm2003-IL2","ChTm2015-IL1","ChTm2015-IL2"))+
  scale_fill_manual(values=fill_mx)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()
####VIOLIN PLOT WR in sympatry with Chalqueño
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/ChTmPeriod_vp.png", width = 600, height = 300)
a <- ggplot2.violinplot(data=Fst_ChTm, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="gray57", meanPointFill="gray57",
                        groupName = 'id',
                        groupColors = fill_mx)+guides(fill=FALSE)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(labels=c("ChTm1984-IL1","ChTm1984-IL2","ChTm2003-IL1",
                            "ChTm2003-IL2","ChTm2015-IL1","ChTm2015-IL2"))
a
dev.off()

### #### Landraces: Zamorano
Fst_ZmNa <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_2_ZmNa", ] #Select the columns with id and Fst value
class(Fst_ZmNa$id)
lr.pops <- levels(as.factor(Fst_ZmNa$id)) #Popnames
lr.pops
#Perform a Mann-Whitney test for Zamorano landrace collected at different periods
Fst_2_ZmNa.stat <- kruskal.test(x=Fst_ZmNa$Fst,g= as.factor(Fst_ZmNa[,1]))
Fst_2_ZmNa.RelabeledSamples.ph <-capture.output(pairwise.wilcox.test(Fst_ZmNa$Fst,
                                                    Fst_ZmNa$id,
                                                    p.adjust.method="bonferroni"),
                                                    file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_2_ZmNa.RelabeledSamples.ph")
####  Boxplot Zamorano
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/ZmNaPeriod_bp.png", width = 700, height = 300)
fill.lr <- c ("dodgerblue4","dodgerblue4","dodgerblue","dodgerblue","blueviolet","blueviolet","firebrick3","firebrick3")
line.lr <- c ("gray57")
Fst_ZmNa_bp <- ggplot(Fst_ZmNa, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)

Fst_ZmNa_bp + scale_x_discrete(labels=c("ZmNa1944-IL1", "ZmNa1944-IL2", 
                                        "ZmNa1979-IL1","ZmNa1979-IL2",
                                        "ZmNa2004-IL1","ZmNa2004-IL2",
                                        "ZmNa2015-IL1","ZmNa2015-IL2"))+
  scale_fill_manual(values=fill.lr)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/ZmNaPeriod_vp.png", width = 700, height = 300)
a <- ggplot2.violinplot(data=Fst_ZmNa, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="black", meanPointFill="black",
                        groupName = 'id',
                        groupColors = fill.lr)+guides(fill=FALSE)+
  theme(axis.text.x  = element_text(size=10),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(labels=c("ZmNa1944-IL1",
                            "ZmNa1944-IL2", 
                            "ZmNa1979-IL1","ZmNa1979-IL2",
                            "ZmNa2004-IL1","ZmNa2004-IL2",
                            "ZmNa2015-IL1","ZmNa2015-IL2"))
a
dev.off()

####  #### WR in sympatry with Zamorano
Fst_ZmTm <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_2_ZmTm", ]
wr.pops <- levels(as.factor(Fst_ZmTm$id)) #Popnames
Fst_2_ZmTm.stat <- kruskal.test(x=Fst_ZmTm$Fst,g= as.factor(Fst_ZmTm[,1]))
capture.output(pairwise.wilcox.test(Fst_ZmTm$Fst,Fst_ZmTm$id,
                                    p.adjust.method="bonferroni"),
               file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_2_ZmTm.RelabeledSamples.ph")
####  Boxplot WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/ZmTmPeriod_bp.png", width = 600, height = 300)
fill_mx <- c("red","red","black","black ","dimgray","dimgray")
line_mx <- "gray57"
Fst_ZmTm_bp <- ggplot(Fst_ZmTm, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)
Fst_ZmTm_bp + scale_x_discrete(labels=c("ZmTm1984-IL1" ,"ZmTm1984-IL2", 
                                        "ZmTm2004-IL1", "ZmTm2004-L2" ,
                                        "ZmTm2015-IL1" ,"ZmTm2015-IL2"))+
  scale_fill_manual(values=fill_mx)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()
####VIOLIN PLOT WR with Zamorano
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/ZmTmPeriod_vp.png", width = 600, height = 300)
a <- ggplot2.violinplot(data=Fst_ZmTm, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="gray57", meanPointFill="gray57",
                        groupName = 'id',
                        groupColors = fill_mx)+guides(fill=FALSE)+
  theme(axis.text.x  = element_text(size=10),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(labels=c("ZmTm1984-IL1" ,"ZmTm1984-IL2", 
                            "ZmTm2004-IL1", "ZmTm2004-IL2" ,
                            "ZmTm2015-IL1" ,"ZmTm2015-IL2"))
a
dev.off()

############################################################################
########  Mushito 
#############################################################################
### #### Landraces: Mushito
Fst_MsNa <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_3_MsNa", ] #Select the columns with id and Fst value
class(Fst_MsNa$id)
lr.pops <- levels(as.factor(Fst_MsNa$id)) #Popnames
lr.pops
#Perform a Mann-Whitney test for Mushito landrace collected at different periods
Fst_3_MsNa.stat <- kruskal.test(x=Fst_MsNa$Fst,g= as.factor(Fst_MsNa[,1]))
Fst_3_MsNa.RelabeledSamples.ph <-capture.output(pairwise.wilcox.test(Fst_MsNa$Fst,
                                                    Fst_MsNa$id,
                                                    p.adjust.method="bonferroni"),
                                                    file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_3_MsNa.RelabeledSamples.ph")
####  Boxplot Mushito
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/MsNaPeriod_bp.png", width = 700, height = 300)
fill.lr <- c ("dodgerblue4","dodgerblue4","dodgerblue","dodgerblue",
              "blueviolet","blueviolet","firebrick3","firebrick3")
line.lr <- c ("gray57")
Fst_MsNa_bp <- ggplot(Fst_MsNa, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)

Fst_MsNa_bp + scale_x_discrete(labels=c("MsNa1952-IL1","MsNa1952-IL2",
                                        "MsNa1979-IL1","MsNa1979-IL2", 
                                        "MsNa2004-IL1","MsNa2004-IL2",
                                        "MsNa2015-IL1","MsNa2015-IL2"))+
  scale_fill_manual(values=fill.lr)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
         size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/MsNaPeriod_vp.png", width = 700, height = 300)
a <- ggplot2.violinplot(data=Fst_MsNa, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="black", meanPointFill="black",
                        groupName = 'id',
                        groupColors = fill.lr)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=11),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(labels=c("MsNa1952-IL1","MsNa1952-IL2",
                            "MsNa1979-IL1","MsNa1979-IL2",
                            "MsNa2004-IL1","MsNa2004-IL2",
                            "MsNa2015-IL1","MsNa2015-IL2"))
a
dev.off()

####  #### WR in sympatry with Mushito
Fst_MsTm <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_3_MsTm", ]
wr.pops <- levels(as.factor(Fst_MsTm$id)) #Popnames
wr.pops
Fst_3_MsTm.stat <- kruskal.test(x=Fst_MsTm$Fst,g= as.factor(Fst_MsTm[,1]))
capture.output(pairwise.wilcox.test(Fst_MsTm$Fst,Fst_MsTm$id,
                                    p.adjust.method="bonferroni"),
               file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_3_MsTm.RelabeledSamples.ph")
####Boxplot WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/MsTmPeriod_bp.png", width = 400, height = 300)
fill_mx <- c("black","black ","dimgray","dimgray")
line_mx <- "gray57"
Fst_MsTm_bp <- ggplot(Fst_MsTm, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)
Fst_MsTm_bp + scale_x_discrete(labels=c("MsTm2004_IL1","MsTm2004_IL2",  
                                        "MsTm2015_IL1","MsTm2015_IL2"))+
  scale_fill_manual(values=fill_mx)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()

png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/MsTmPeriod_vp.png", width = 400, height = 300)
Fst_MsTm_vp<- ggplot2.violinplot(data=Fst_MsTm, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="gray57", meanPointFill="gray57",
                        groupName = 'id',
                        groupColors = fill_mx)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ 
  scale_x_discrete(labels=c("MsTm2004_IL1","MsTm2004_IL2",  
                            "MsTm2015_IL1","MsTm2015_IL2"))
Fst_MsTm_vp
dev.off()
### #### Landraces: Tabloncillo
Fst_TbNa <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_4_TbNa", ] #Select the columns with id and Fst value
class(Fst_TbNa$id)
lr.pops <- levels(as.factor(Fst_TbNa$id)) #Popnames
lr.pops
#Perform a Mann-Whitney test for Tabloncillo landrace collected at different periods
Fst_4_TbNa.stat <- kruskal.test(x=Fst_TbNa$Fst,g= as.factor(Fst_TbNa[,1]))
Fst_4_TbNa.RelabeledSamples.ph <-capture.output(pairwise.wilcox.test(Fst_TbNa$Fst,
                                                   Fst_TbNa$id,
                                                   p.adjust.method="bonferroni"),
                                                file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_4_TbNa.RelabeledSamples.ph")

####  Boxplot Tabloncillo
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/TbNaPeriod_bp.png", width = 900, height = 300)
fill.lr <- c ("dodgerblue4","dodgerblue4","dodgerblue","dodgerblue","cyan","cyan","blueviolet","blueviolet","firebrick3","firebrick3")
line.lr <- c ("gray57")
Fst_TbNa_bp <- ggplot(Fst_TbNa, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)
Fst_TbNa_bp + scale_x_discrete(labels=c("TbNa1946a-IL1","TbNa1946a-IL2",
                                        "TbNa1946b-IL1","TbNa1946b-IL2",
                                        "TbNa1960-IL1","TbNa1960-IL2",
                                        "TbNa2003-IL1","TbNa2003-IL2",
                                        "TbNa2015-IL1","TbNa2015-IL2"))+
  scale_fill_manual(values=fill.lr)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()


png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/TbNaPeriod_vp.png", width = 900, height = 300)
a <- ggplot2.violinplot(data=Fst_TbNa, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="black", meanPointFill="black",
                        groupName = 'id',
                        groupColors = fill.lr)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ 
  scale_x_discrete(labels=c("TbNa1946a-IL1","TbNa1946a-IL2" ,
                            "TbNa1946b_IL1","TbNa1946b-IL2",
                            "TbNa1960-IL1","TbN1960-IL2",
                            "TbNa2003-IL1","TbNa2003-IL2",
                            "TbNa2015-IL1","TbNa2015-IL2"))
a
dev.off()

####  #### WR in sympatry with Tabloncillo
Fst_TbTp<- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_4_TbTp", ]
wr.pops <- levels(as.factor(Fst_TbTp$id)) #Popnames
Fst_4_TbTp.stat <- kruskal.test(x=Fst_TbTp$Fst,g= as.factor(Fst_TbTp[,1]))
capture.output(pairwise.wilcox.test(Fst_TbTp$Fst,Fst_TbTp$id,
                                    p.adjust.method="bonferroni"),file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_4_TbTp.RelabeledSamples.ph")
####Boxplot WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/TbTpPeriod_bp.png", width = 600, height = 300)
fill_mx <- c("red","red","black","black ","dimgray","dimgray")
line_mx <- "gray57"
Fst_TbTp_bp <- ggplot(Fst_TbTp, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)
Fst_TbTp_bp + scale_x_discrete(labels=c ("TbTp1983-IL1","TbTp1983-IL2",
                                         "TbTp2003-IL1","TbTp2003-IL2",
                                         "TbTp2015-IL1","TbTp2015-IL2"))+
  scale_fill_manual(values=fill_mx)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()
####VIOLIN PLOT WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/TbTpPeriod_vp.png", width = 600, height = 300)
a <- ggplot2.violinplot(data=Fst_TbTp, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="gray57", meanPointFill="gray57",
                        groupName = 'id',
                        groupColors = fill_mx)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ 
  scale_x_discrete(labels=c("TbTp1983-IL1","TbTp1983-IL2",
                            "TbTp2003-IL1","TbTp2003-IL2",
                            "TbTp2015-IL1","TbTp2015-IL2"))
a
dev.off()





### #### Landraces: Pepitilla
#################################
Fst_PpNa <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_5_PpNa", ] #Select the columns with id and Fst value
class(Fst_PpNa$id)
lr.pops <- levels(as.factor(Fst_PpNa$id)) #Popnames
lr.pops
#Perform a Mann-Whitney test for Pepitilla landrace collected at different periods
Fst_5_PpNa.stat <- kruskal.test(x=Fst_PpNa$Fst,g= as.factor(Fst_PpNa[,1]))
Fst_5_PpNa.RelabeledSamples.ph <-capture.output(pairwise.wilcox.test(Fst_PpNa$Fst,
                                                Fst_PpNa$id,
                                                p.adjust.method="bonferroni"),
                                                file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_5_PpNa.RelabeledSamples.ph")
####  Boxplot Pepitilla
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/PpNaPeriod_bp.png", width = 900, height = 300)
fill.lr <- c ("dodgerblue4","dodgerblue4","dodgerblue","dodgerblue","cyan","cyan","blueviolet","blueviolet","firebrick3","firebrick3")
line.lr <- c ("gray57")
Fst_PpNa_bp <- ggplot(Fst_PpNa, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)

Fst_PpNa_bp + scale_x_discrete(labels=c("PpNa1943-IL1","PpNa1943-IL2",
                                        "PpNa1960-IL1","PpNa1960-IL2",
                                        "PpNa1973-IL1","PpNa1973-IL2",
                                        "PpNa2003a-IL1","PpNa2003a-IL2",
                                        "PpNa2003b-IL1","PpNa2003b-IL2"))+
  scale_fill_manual(values=fill.lr)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/PpNaPeriod_vp.png", width = 900, height = 300)
a <- ggplot2.violinplot(data=Fst_PpNa, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="black", meanPointFill="black",
                        groupName = 'id',
                        groupColors = fill.lr)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+ 
  scale_x_discrete(labels=c("PpNa1943-IL1","PpNa1943-IL2",
                            "PpNa1960-IL1","PpNa1960-IL2",
                            "PpNa1973-IL1","PpNa1973-IL2",
                            "PpNa2003a-IL1","PpNa2003a-IL2",
                            "PpNa2003b-IL1","PpNa2003b-IL2"))
a
dev.off()

####  #### WR in sympatry with Pepitilla
Fst_PpTp<- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_5_PpTp", ]
wr.pops<- levels(as.factor(Fst_PpTp$id)) #Popnames
wr.pops
Fst_5_PpTp.stat <- kruskal.test(x=Fst_PpTp$Fst,g= as.factor(Fst_PpTp[,1]))
capture.output(pairwise.wilcox.test(Fst_PpTp$Fst,Fst_PpTp$id,
                                    p.adjust.method="bonferroni"),
                                    file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_5_PpTp.RelabeledSamples.ph")
####Boxplot WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/PpTpPeriod_bp.png", width = 700, height = 300)
fill_mx <- c("red","red","black","black ","dimgray","dimgray","bisque3","bisque3")
line_mx <- "gray57"
Fst_PpTp_bp <- ggplot(Fst_PpTp, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)
Fst_PpTp_bp + scale_x_discrete(labels=c ("PpTp1984-IL1","PpTp1984-IL2",
                                         "PpTp2003a-IL1","PpTp2003a-IL2", 
                                         "PpTp2003b-IL1","PpTp2003b-IL2",
                                         "PpTp2011-IL1","PpTp2011-IL2"))+
  scale_fill_manual(values=fill_mx)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()
####VIOLIN PLOT WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/PpTpPeriod_vp.png", width = 700, height = 300)
a <- ggplot2.violinplot(data=Fst_PpTp, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="gray57", meanPointFill="gray57",
                        groupName = 'id',
                        groupColors = fill_mx)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_discrete(labels=c("PpTp1984-IL1","PpTp1984-IL2",
                            "PpTp2003a-IL1","PpTp2003a-IL2", 
                            "PpTp2003b-IL1","PpTp2003b-IL2",
                            "PpTp2011-IL1","PpTp2011-IL2"))
a
dev.off()



### #### Landraces: Vandeño
#############################
Fst_VaNa <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_6_VaNa", ] #Select the columns with id and Fst value
class(Fst_VaNa$id)
lr.pops <- levels(as.factor(Fst_VaNa$id)) #Popnames
lr.pops
#Perform a Mann-Whitney test for Vandeño landrace collected at different periods
Fst_6_VaNa.stat <- kruskal.test(x=Fst_VaNa$Fst,g= as.factor(Fst_VaNa[,1]))
Fst_6_VaNa.RelabeledSamples.ph <-capture.output(pairwise.wilcox.test(Fst_VaNa$Fst,
                                                    Fst_VaNa$id,
                                                    p.adjust.method="bonferroni"),file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_6_VaNa.RelabeledSamples.ph")
####  Boxplot Vandeño
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/VaNaPeriod_bp.png", width = 600, height = 300)
fill.lr <- c ("dodgerblue4","dodgerblue4","blueviolet","blueviolet","firebrick3","firebrick3")
line.lr <- c ("gray57")
Fst_VaNa_bp <- ggplot(Fst_VaNa, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)

Fst_VaNa_bp + scale_x_discrete(labels=c("VaNa1943-IL1","VaNa1943-IL2",
                                        "VaNa2003-IL1","VaNa2003-IL2",
                                        "VaNa2015-IL1","VaNa2015-IL2"))+
  scale_fill_manual(values=fill.lr)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/VaNaPeriod_vp.png", width = 600, height = 300)
a <- ggplot2.violinplot(data=Fst_VaNa, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="black", meanPointFill="black",
                        groupName = 'id',
                        groupColors = fill.lr)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_discrete(labels=c("VaNa1943-IL1","VaNa1943-IL2",
                            "VaNa2003-IL1","VaNa2003-IL2",
                            "VaNa2015-IL1","VaNa2015-IL2"))
a
dev.off()

####  #### WR in sympatry with Vandeño
Fst_VaTp<- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_6_VaTp", ]
wr.pops<- levels(as.factor(Fst_VaTp$id)) #Popnames
wr.pops
Fst_6_VaTp.stat <- kruskal.test(x=Fst_VaTp$Fst,g= as.factor(Fst_VaTp[,1]))
capture.output(pairwise.wilcox.test(Fst_VaTp$Fst,Fst_VaTp$id,
                                    p.adjust.method="bonferroni"),file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_6_VaTp.RelabeledSamples.ph")
####Boxplot WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/VaTpPeriod_bp.png", width = 200, height = 300)
fill_mx <- c("dimgray","dimgray","bisque3","bisque3")
line_mx <- "gray57"
Fst_VaTp_bp <- ggplot(Fst_VaTp, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)
Fst_VaTp_bp + scale_x_discrete(labels=c ("VaTp2003-IL1" ,"VaTp2003-IL2"))+
  scale_fill_manual(values=fill_mx)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()
####VIOLIN PLOT WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/VaTpPeriod_vp.png", width = 200, height = 300)
a <- ggplot2.violinplot(data=Fst_VaTp, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="gray57", meanPointFill="gray57",
                        groupName = 'id',
                        groupColors = fill_mx)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_discrete(labels=c("VaTp2003-IL1" ,"VaTp2003-IL2"))
a
dev.off()

###################################################################################################################
### #### Landraces: Bolita
Fst_BoNa <- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_7_BoNa", ] #Select the columns with id and Fst value
class(Fst_BoNa$id)
lr.pops <- levels(as.factor(Fst_BoNa$id)) #Popnames
lr.pops
#Perform a Mann-Whitney test for Bolita landrace collected at different periods
Fst_6_BoNa.stat <- kruskal.test(x=Fst_BoNa$Fst,g= as.factor(Fst_BoNa[,1]))
Fst_6_BoNa.RelabeledSamples.ph <-capture.output(pairwise.wilcox.test(Fst_BoNa$Fst,
                                                Fst_BoNa$id,
                                                p.adjust.method="bonferroni"),
                                                file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_7_BoNa.RelabeledSamples.ph")
####  Boxplot Bolita
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/BoNaPeriod_bp.png", width = 900, height = 300)
fill.lr <- c ("dodgerblue4","dodgerblue4","dodgerblue","dodgerblue",
              "cyan","cyan","blueviolet","blueviolet","firebrick3","firebrick3")
line.lr <- c ("gray57")
Fst_BoNa_bp <- ggplot(Fst_BoNa, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)

Fst_BoNa_bp + scale_x_discrete(labels=c( "BoNa1943-IL1","BoNa1943-IL2",
                                         "BoNa1970-IL1","BoNa1970-IL2",  
                                         "BoNa2003-IL1","BoNa2003-IL2",
                                         "BoNa2015a-IL1","BoNa2015a-IL2", 
                                         "BoNa2015b-IL1","BoNa2015b-IL2"))+
  scale_fill_manual(values=fill.lr)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
        size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/BoNaPeriod_vp.png", width = 900, height = 300)
a <- ggplot2.violinplot(data=Fst_BoNa, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="black", meanPointFill="black",
                        groupName = 'id',
                        groupColors = fill.lr)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_discrete(labels=c("BoNa1943-IL1","BoNa1943-IL2",
                            "BoNa1970-IL1","BoNa1970-IL2",  
                            "BoNa2003-IL1","BoNa2003-IL2",
                            "BoNa2015a-IL1","BoNa2015a-IL2", 
                            "BoNa2015b-IL1","BoNa2015b-IL2"))
a
dev.off()

####  #### WR in sympatry with Bolita
Fst_BoTp<- Fst_5kb[,c(1,8)][Fst_5kb$id %like% "Fst_7_BoTp", ]
wr.pops<- levels(as.factor(Fst_BoTp$id)) #Popnames
wr.pops
Fst_7_BoTp.stat <- kruskal.test(x=Fst_BoTp$Fst,g= as.factor(Fst_BoTp[,1]))
capture.output(pairwise.wilcox.test(Fst_BoTp$Fst,Fst_BoTp$id,
                                    p.adjust.method="bonferroni"),file = "~/hubic/MisScripts/15.Fst_Windows/out/Stats/Fst_7_BoTp.RelabeledSamples.ph")
####Boxplot WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/BoTpPeriod_bp.png", width = 400, height = 300)
fill_mx <- c("dimgray","dimgray","bisque3","bisque3")
line_mx <- "gray57"
Fst_BoTp_bp <- ggplot(Fst_BoTp, aes(x = id, y = Fst, fill= id)) +
  geom_boxplot(aes(fill=id))+ guides(fill=FALSE)
Fst_BoTp_bp + scale_x_discrete(labels=c ("BoTp2003_IL1","BoTp2003_IL2",
                                         "BoTp2015_IL1","BoTp2015_IL2"))+
  scale_fill_manual(values=fill_mx)+
  theme(text = element_text(size=13),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
        size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()
####VIOLIN PLOT WR
png("~/hubic/MisScripts/15.Fst_Windows/out/Plots/RelabeledSamples/BoTpPeriod_vp.png", width = 400, height = 300)
a <- ggplot2.violinplot(data=Fst_BoTp, xName='id',yName='Fst',
                        addMean=TRUE,meanPointShape=20, meanPointSize=3,
                        meanPointColor="gray57", meanPointFill="gray57",
                        groupName = 'id',
                        groupColors = fill_mx)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size=9),
        axis.title.x = element_blank(),
        axis.line = element_line(colour = "black", 
                                 size = 0.5, linetype = "solid"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_discrete(labels=c("BoTp2003_IL1","BoTp2003_IL2",
                            "BoTp2015_IL1","BoTp2015_IL2"))
a
dev.off()

