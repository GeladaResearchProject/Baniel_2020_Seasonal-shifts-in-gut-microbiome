rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(phyloseq)
library(tidyverse)
library(vegan)
library(plyr) 
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
theme_set (theme_classic(base_size=20))



############## DATASET #################
load("~/Desktop/gelada_SILVA_raref_ASV filter.RData")
gelada_physeq #758 samples - 32853 ASVs
otu_table(gelada_physeq)[1:5,1:5] #rows=ASVs,columns=samples
head(sample_data(gelada_physeq))
summary(sample_data(gelada_physeq)$NbReads) #20K reads randomly selected in the fastq files then fitlered
meta<-data.frame(sample_data(gelada_physeq))




################################### 1. ALPHA DIVERSITY ####################################
### Calculate alpha diversity indices ###
alpha<-estimate_richness(gelada_physeq,measures = c("Observed","Shannon"))

#Faith PD
library(picante)
faith<-pd(as.matrix(as.data.frame(t(otu_table(gelada_physeq)))),phy_tree(gelada_physeq),include.root=TRUE)#row=samples,columns=ASV
rownames(faith)==rownames(alpha)
alpha1<-cbind(rownames(alpha),alpha,faith)
colnames(alpha1)[1]<-"SampleID"
head(alpha1)

#Add metadata
alpha2<-merge(alpha1,meta,by="SampleID",all.x=TRUE,all.y=FALSE)
head(alpha2)


### Mixed models ###

#z-transformation of quantitative variables
stdize=function(x) {(x - mean(x))/sd(x)}
alpha2$zAgeYear<-stdize(alpha2$AgeYear)
alpha2$zNbReads<-stdize(alpha2$NbReads)
alpha2$zRain30<-stdize(alpha2$Rain30)
alpha2$zMinT30<-stdize(alpha2$MinT30)

#LMM
m1<-lmer(Shannon ~ Sex + zAgeYear + zRain30 + zMinT30 + zNbReads + (1|Unit) + (1|ID), data=alpha2, na.action=na.fail)
summary(m1)
drop1(m1,test="Chisq")
confint(m1,level=0.95,method="boot")
par(mfrow=c(2,2),mar=c(5,5,5,5))
visreg(m1,"Sex",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,cex.main=2,ylab="Observed richness",main="Partial residuals")
visreg(m1,"zAgeYear",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Observed richness")
visreg(m1,"zRain30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Observed richness")
visreg(m1,"zMinT30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Observed richness")

m2<-lmer(Observed ~ Sex + zAgeYear + zRain30 + zMinT30 + zNbReads + (1|Unit) + (1|ID), data=alpha2, na.action=na.fail)
summary(m2)
drop1(m2,test="Chisq")
confint(m2,level=0.95,method="boot")
visreg(m2,"Sex",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,cex.main=2,ylab="Shannon index",main="Partial residuals")
visreg(m2,"zAgeYear",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Shannon index")
visreg(m2,"zRain30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Shannon index")
visreg(m2,"zMinT30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Shannon index")


m3<-lmer(PD ~ Sex + zAgeYear + zRain30 + zMinT30 + zNbReads + (1|Unit) + (1|ID), data=alpha2, na.action=na.fail)
summary(m3)
drop1(m3,test="Chisq")
confint(m3,level=0.95,method="boot")
visreg(m3,"Sex",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,cex.main=2,ylab="Faith's PD",main="Partial residuals")
visreg(m3,"zAgeYear",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Faith's PD")
visreg(m3,"zRain30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Faith's PD")
visreg(m3,"zMinT30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Faith's PD")





########################################## 2. BETA DIVERSITY ##############################################
# Calculate Bray-Curtis distance
bcdist = phyloseq::distance(gelada_physeq, method="bray",normalized=TRUE, parallel=FALSE, fast=TRUE)  
PC1 <- ordinate(gelada_physeq, method = "MDS", distance = bcdist)
meta$bc.PC1=PC1$vectors[,1]
meta$bc.PC2=PC1$vectors[,2]


# Calculate unweighted UniFrac distance
udist=UniFrac(gelada_physeq,weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)
PC2 <- ordinate(gelada_physeq, method = "MDS", distance = udist)
meta$uu.PC1=PC2$vectors[,1]
meta$uu.PC2=PC2$vectors[,2]


# Calculate weighted UniFrac distance
wdist=UniFrac(gelada_physeq, weighted=TRUE, normalized=TRUE, parallel=FALSE, fast=TRUE)
pcoa3 <- ordinate(gelada_physeq, method = "MDS", distance = wdist)
meta$wu.PC1=PC2$vectors[,1]
meta$wu.PC2=PC2$vectors[,2]


#RAIN30
a<-ggplot(meta,aes(x=bc.PC1, y=Rain30))+ ggtitle("Bray-Curtis")+ geom_point() + coord_flip() + ylab("Cumulative rainfall") + xlab("PC1 axis") + theme(legend.position="none")
b<-ggplot(meta,aes(x=bc.PC2, y=Rain30))+ geom_point() + coord_flip() + ylab("Cumulative rainfall") + xlab("PC2 axis") + theme(legend.position="none")
c<-ggplot(meta,aes(x=uu.PC1, y=Rain30))+ ggtitle("Unweighted UniFrac")+ geom_point() + coord_flip() + ylab("Cumulative rainfall") + xlab("PC1 axis") + theme(legend.position="none")
d<-ggplot(meta,aes(x=uu.PC2, y=Rain30))+ geom_point() + coord_flip() + ylab("Cumulative rainfall") + xlab("PC2 axis") + theme(legend.position="none")
e<-ggplot(meta,aes(x=wu.PC1, y=Rain30))+ ggtitle("Weighted UniFrac")+ geom_point() + coord_flip() + ylab("Cumulative rainfall") + xlab("PC1 axis") + theme(legend.position="none")
f<-ggplot(meta,aes(x=wu.PC2, y=Rain30))+ geom_point() + coord_flip() + ylab("Cumulative rainfall") + xlab("PC2 axis") + theme(legend.position="none")
grid.arrange(nrow=3,ncol=2,a,b,c,d,e,f)


#MinT30
a<-ggplot(meta,aes(x=bc.PC1, y=MinT30))+ ggtitle("Bray-Curtis")+ geom_point() + coord_flip() + ylab("Minimum temeprature") + xlab("PC1 axis") + theme(legend.position="none")
b<-ggplot(meta,aes(x=bc.PC2, y=MinT30))+ geom_point() + coord_flip() + ylab("Minimum temeprature") + xlab("PC2 axis") + theme(legend.position="none")
c<-ggplot(meta,aes(x=uu.PC1, y=MinT30))+ ggtitle("Unweighted UniFrac")+ geom_point() + coord_flip() + ylab("Minimum temeprature") + xlab("PC1 axis") + theme(legend.position="none")
d<-ggplot(meta,aes(x=uu.PC2, y=MinT30))+ geom_point() + coord_flip() + ylab("Minimum temeprature") + xlab("PC2 axis") + theme(legend.position="none")
e<-ggplot(meta,aes(x=wu.PC1, y=MinT30))+ ggtitle("Weighted UniFrac")+ geom_point() + coord_flip() + ylab("Minimum temeprature") + xlab("PC1 axis") + theme(legend.position="none")
f<-ggplot(meta,aes(x=wu.PC2, y=MinT30))+ geom_point() + coord_flip() + ylab("Minimum temeprature") + xlab("PC2 axis") + theme(legend.position="none")
grid.arrange(nrow=3,ncol=2,a,b,c,d,e,f)

#AGE
a<-ggplot(meta,aes(x=bc.PC1, y=AgeYear))+ ggtitle("Bray-Curtis")+ geom_point() + coord_flip() + ylab("Age (in years)") + xlab("PC1 axis") + theme(legend.position="none")
b<-ggplot(meta,aes(x=bc.PC2, y=AgeYear))+ geom_point() + coord_flip() + ylab("Age (in years)") + xlab("PC2 axis") + theme(legend.position="none")
c<-ggplot(meta,aes(x=uu.PC1, y=AgeYear))+ ggtitle("Unweighted UniFrac")+ geom_point() + coord_flip() + ylab("Age (in years)") + xlab("PC1 axis") + theme(legend.position="none")
d<-ggplot(meta,aes(x=uu.PC2, y=AgeYear))+ geom_point() + coord_flip() + ylab("Age (in years)") + xlab("PC2 axis") + theme(legend.position="none")
e<-ggplot(meta,aes(x=wu.PC1, y=AgeYear))+ ggtitle("Weighted UniFrac")+ geom_point() + coord_flip() + ylab("Age (in years)") + xlab("PC1 axis") + theme(legend.position="none")
f<-ggplot(meta,aes(x=wu.PC2, y=AgeYear))+ geom_point() + coord_flip() + ylab("Age (in years)") + xlab("PC2 axis") + theme(legend.position="none")
grid.arrange(nrow=3,ncol=2,a,b,c,d,e,f)


#SEX
a<-ggplot(meta,aes(x=bc.PC1,y=bc.PC2, colour=Sex))+ ggtitle("Bray-Curtis")+ geom_point() + xlab("PC1 axis") + ylab("PC2 axis") 
b<-ggplot(meta,aes(x=uu.PC1,y=uu.PC2, colour=Sex))+ ggtitle("Unweighted UniFrac")+ geom_point() + xlab("PC1 axis") + ylab("PC2 axis") 
c<-ggplot(meta,aes(x=wu.PC1,y=wu.PC2, colour=Sex))+ ggtitle("Weighted UniFrac")+ geom_point() + xlab("PC1 axis") + ylab("PC2 axis") 
grid.arrange(nrow=1,a,b,c)




### PERMANOVA test ####
#Effect of individuals
i1<-adonis2(bcdist ~ ID + NbReads, data = meta, permutations = 10000)
i2<-adonis2(udist ~ ID + NbReads, data = meta, permutations = 10000)
i3<-adonis2(wdist ~ ID + NbReads, data = meta, permutations = 10000)

#Bray Curtis
adonis2(bcdist ~ NbReads + Unit +  Rain30 + MinT30 + Sex + AgeYear, strata=ID, data = meta, permutations = 10000)

#Unweighted UniFrac
adonis2(udist ~ NbReads + Unit +  Rain30 + MinT30 + Sex + AgeYear, strata=ID, data = meta, permutations = 10000)

#Weighted UniFrac
adonis2(wdist ~NbReads + Unit +  Rain30 + MinT30 + Sex + AgeYear, strata=ID, data = meta, permutations = 10000)


