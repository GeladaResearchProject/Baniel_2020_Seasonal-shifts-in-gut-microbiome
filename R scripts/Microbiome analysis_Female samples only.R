rm(list=ls())
'%!in%'<-function(x,y)!('%in%'(x,y))
library(phyloseq)
library(tidyverse)
library(vegan)
library(plyr) 
library(reshape2)
library(lme4)
library(pbkrtest)
library(DHARMa)
library(visreg)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(multcomp)
theme_set (theme_classic(base_size=20))



############## DATASET #################
load("~/Desktop/gelada_SILVA_ASV filter.RData")
gelada_physeq #758 samples - 3295 ASVs
otu_table(gelada_physeq)[1:5,1:5]#rows=ASVs,columns=samples
head(sample_data(gelada_physeq))

# Restrict to female samples
gelada_physeq = prune_samples(sample_data(gelada_physeq)$ReproState %in% c("C","P","L1"), gelada_physeq)
gelada_physeq #439 samples

meta<-data.frame(sample_data(gelada_physeq))
head(meta)

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
m1<-lmer(Shannon ~ ReproState + zAgeYear + zRain30 + zMinT30 + zNbReads + (1|Unit) + (1|ID), data=alpha2, na.action=na.fail)
summary(m1)
drop1(m1,test="Chisq")
confint(m1,level=0.95,method="boot")
par(mfrow=c(2,2),mar=c(5,5,5,5))
visreg(m1,"ReproState",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,cex.main=2,ylab="Observed richness",main="Partial residuals")
visreg(m1,"zAgeYear",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Observed richness")
visreg(m1,"zRain30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Observed richness")
visreg(m1,"zMinT30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Observed richness")

m2<-lmer(Observed ~ ReproState + zAgeYear + zRain30 + zMinT30 + zNbReads + (1|Unit) + (1|ID), data=alpha2, na.action=na.fail)
summary(m2)
drop1(m2,test="Chisq")
confint(m2,level=0.95,method="boot")
visreg(m2,"ReproState",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,cex.main=2,ylab="Shannon index",main="Partial residuals")
visreg(m2,"zAgeYear",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Shannon index")
visreg(m2,"zRain30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Shannon index")
visreg(m2,"zMinT30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Shannon index")


m3<-lmer(PD ~ ReproState + zAgeYear + zRain30 + zMinT30 + zNbReads + (1|Unit) + (1|ID), data=alpha2, na.action=na.fail)
summary(m3)
drop1(m3,test="Chisq")
confint(m3,level=0.95,method="boot")
visreg(m3,"ReproState",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,cex.main=2,ylab="Faith's PD",main="Partial residuals")
visreg(m3,"zAgeYear",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Faith's PD")
visreg(m3,"zRain30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Faith's PD")
visreg(m3,"zMinT30",font.lab=2,points=list(cex=1, pch=16,col="black"),cex.lab=1.5,cex.axis=1,ylab="Faith's PD")



########################################## 2. BETA DIVERSITY ##############################################
m1<- as.matrix(t(otu_table(gelada_physeq)))
m1[1:5,1:5]#row=samples
taxa<-data.frame(tax_table(gelada_physeq))
taxa$taxa_id<-rownames(taxa)

### clr-transfornmation of the counts ###
#for more information on compositional microbiome analysis, see Gloor et al. 2017 and Gloor's github tutorial (https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot)
pseudocount <- 0.65 #add a pseudo count for zeroes
clr <- t(apply(t(m1)+pseudocount, 2, compositions::clr))
clr[1:5,1:5] #rows=samples
dist(clr) #Aitchison distance

### Perform the PCA ###
pca <- prcomp(clr)
plot(pca$x[,1],pca$x[,2])

# Calculate the variance explained by PC1 and PC2
d.mvar <- sum(pca$sdev^2) # total variance
PC1 <- paste("PC1 ","(",round(sum(pca$sdev[1]^2)/d.mvar, 3)*100,"%",")",sep="")
PC2 <- paste("PC2 ","(",round(sum(pca$sdev[2]^2)/d.mvar, 3)*100,"%",")",sep="")

#Add metadata
df<-data.frame(pca$x)
df$SampleID<-rownames(df)
df1<-merge(df,meta,by="SampleID",all.x=TRUE,all.y=FALSE)
names(df1)
df1[1:5,c(1:3,441:454)]

#PCA representation
a<-ggplot(df1, aes(x=PC1,y=PC2,colour=ReproState)) + geom_point()+ xlab(PC1)+ ylab(PC2)+theme(axis.title=element_text(size=18,face="bold"))
b<-ggplot(df1, aes(x=Sex,y=PC1,fill=ReproState)) + geom_boxplot()+ ylab(PC1)+xlab("Reproductive State")+theme(axis.title=element_text(size=18,face="bold"))+theme(legend.position="none")
c<-ggplot(df1, aes(x=Sex,y=PC2,fill=ReproState)) + geom_boxplot()+ ylab(PC2)+xlab("Reproductive State")+theme(axis.title=element_text(size=18,face="bold"))+theme(legend.position="none")
ggarrange(a,b,c,nrow=1)

a<-ggplot(df1, aes(x=Rain30,y=PC1)) + geom_point() + ylab(PC1) + xlab("Cumulative rainfall") + theme(axis.title=element_text(size=18,face="bold"))+scale_x_continuous(breaks=seq(0,600,100),limits=c(0,600))
b<-ggplot(df1, aes(x=Rain30,y=PC2)) + geom_point() + ylab(PC1) + xlab("Cumulative rainfall") + theme(axis.title=element_text(size=18,face="bold"))+scale_x_continuous(breaks=seq(0,600,100),limits=c(0,600))
ggarrange(a,b,nrow=1)

a<-ggplot(df1, aes(x=MinT30,y=PC1)) + geom_point() + ylab(PC1) + xlab("Minimum temperature") + theme(axis.title=element_text(size=18,face="bold"))+scale_x_continuous(breaks=seq(5.8,12,2),limits=c(5.8,12))
b<-ggplot(df1, aes(x=MinT30,y=PC2)) + geom_point() + ylab(PC2) + xlab("Minimum temperature") + theme(axis.title=element_text(size=18,face="bold"))+scale_x_continuous(breaks=seq(5.8,12,2),limits=c(5.8,12))
ggarrange(a,b,nrow=1)

a<-ggplot(df1, aes(x=AgeYear,y=PC1)) + geom_point()+ ylab(PC1)+xlab("Age (in years)")+theme(axis.title=element_text(size=18,face="bold")) + scale_x_continuous(breaks=seq(5,30,5),limits=c(4.5,30.2))
b<-ggplot(df1, aes(x=AgeYear,y=PC2)) + geom_point()+ ylab(PC2)+xlab("Age (in years)")+theme(axis.title=element_text(size=18,face="bold")) + scale_x_continuous(breaks=seq(5,30,5),limits=c(4.5,30.2))
ggarrange(a,b,nrow=1)


### PERMANOVA test ####
aich<-dist(clr) #Aitchison distance

#Effect of individuals
adonis2(aich ~ ID + NbReads, data = meta, permutations = 10000)

#Effects of predictors
adonis2(aich ~ NbReads + Unit + Rain30 + MinT30 + ReproState + AgeYear, strata=ID, data = meta, permutations = 10000)



###################################### 3. DIFFERENTIAL ABONDANCE TESTING  ################################################
### Clean unknown ASV names from SILVA taxonomy ###
table(tax_table(gelada_physeq)[,"Phylum"],useNA="always")
table(tax_table(gelada_physeq)[,"Class"],useNA="always")
table(tax_table(gelada_physeq)[,"Order"],useNA="always")
table(tax_table(gelada_physeq)[,"Family"],useNA="always")
table(tax_table(gelada_physeq)[,"Genus"],useNA="always")
table(tax_table(gelada_physeq)[,"Species"],useNA="always")
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("Unassigned*",NA,x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("uncultured*",NA, x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("unidentified*",NA, x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("gut metagenome",NA, x)))
tax_table(gelada_physeq)<-apply(tax_table(gelada_physeq),2, function(x) (gsub("metagenome",NA, x)))
tax_table(gelada_physeq)[,"Family"] [tax_table(gelada_physeq)[,"Order"] %in% "WCHB1-41" ]<- "RFP12" 


### Agglomeration ###

# CHOOSE WANTED LEVEL HERE #
ps1 = tax_glom(gelada_physeq,"Phylum",NArm=FALSE)
ps2<-psmelt(ps1)
ps2$Taxa<-ps2$Phylum

ps1 = tax_glom(gelada_physeq,"Class",NArm=FALSE)
ps2<-psmelt(ps1)
ps2$Taxa<-ps2$Class

ps1 = tax_glom(gelada_physeq,"Order",NArm=FALSE)
ps2<-psmelt(ps1)
ps2$Taxa<-ps2$Order

ps1 = tax_glom(gelada_physeq,"Family",NArm=FALSE)
ps2<-psmelt(ps1)
ps2$Taxa<-ps2$Family

ps1 = tax_glom(gelada_physeq,"Genus",NArm=FALSE)
ps2<-psmelt(ps1)
ps2$Taxa<-ps2$Genus


#Calculate relative abundance of each taxa across all samples
head(ps2)
dim(ps2)
ps2$RelativeAbundance<-(as.numeric(as.character(ps2$Abundance))*100)/as.numeric(as.character(ps2$NbReads))
a1<-sapply(unique(ps2$Taxa),function(x) mean(ps2$RelativeAbundance[ps2$Taxa==x],na.rm=TRUE))
a2<-data.frame(cbind(Taxa=as.character(unique(ps2$Taxa)),RA=a1))
a2$RA<-as.numeric(as.character(a2$RA))
a2[with(a2, order(-RA)), ] 


#Filter taxa>0.01% relative abundance for the models
ps3<-droplevels(subset(ps2,ps2$Taxa %!in% NA))
ps3<-droplevels(subset(ps3,ps3$Taxa %in% unique(a2$Taxa[a2$RA>=0.01])))
length(unique(ps3$Taxa))
length(unique(ps2$Taxa))


### MODELS ###

#Z-transformation of numerical variables
stdize=function(x) {(x - mean(x))/sd(x)}
ps3$zAgeYear<-stdize(ps3$AgeYear)
ps3$zRain30<-stdize(ps3$Rain30)
ps3$zMinT30<-stdize(ps3$MinT30)
ps3$zNbReads<-stdize(ps3$NbReads)

ps4<-ps3 

#Remove taxa that do not converge with a NB model
#Phylum+Class: all phyla converge 

#Order
ps4<-droplevels(subset(ps4,ps4$Taxa %!in% c("Verrucomicrobiales","Enterobacteriales","Bacillales")))

#Family
ps4<-droplevels(subset(ps4,ps4$Taxa %!in% c("Bacteroidaceae","Akkermansiaceae","Enterobacteriaceae","Staphylococcaceae","Marinifilaceae")))

#Genus
ps4<-droplevels(subset(ps4,ps4$Taxa %!in% c("Oxalobacter","Alistipes","Butyricimonas","Bacteroides","Akkermansia","[Ruminococcus] gnavus group","Adlercreutzia","Escherichia-Shigella","Staphylococcus","Hungatella")))


### LOOP ###
#Fit a NB GLMM for each taxa
df2<-NULL
for (i in unique(ps4$Taxa)){
 data<-droplevels(subset(ps4,ps4$Taxa %in% i))
  m1 <- glmer.nb (Abundance ~ ReproState + zAgeYear + zMinT30 + zRain30 + offset(log(NbReads)) + (1|ID) + (1|Unit),data=data)
  shap<-round(as.numeric(testUniformity(simulateResiduals(m1),plot=F)[2]),4)
  var1<-as.data.frame(print(VarCorr(m1),comp=c("Variance","Std.Dev.")))
  varRE<-round(as.numeric(var1[,4]),3)
  sdRE<-round(as.numeric(var1[,5]),3)
  test = glht(m1,linfct=mcp(ReproState="Tukey"))
  summary(test)
  pval1<-round(as.numeric(summary(test,test=adjusted(type="none"))$test$pvalues),4)
  coef1<-round(as.numeric(summary(test,test=adjusted(type="none"))$test$coef),4)
  coef<-round(summary(m1)$coef[4:6,1],4)
  drop<-drop1(m1,test="Chisq")
  pval<-round(drop[3:5,4],4)
  df<-c(i,shap,varRE,sdRE,coef1,pval1,coef,pval)
  df2<-rbind(df2,df)
  colnames(df2)<-c("Taxa","pShap","varRE.ID","varRE.Unit","sdRE.ID","sdRE.Unit","coefLC","coefPC","coefPL","pLC","pPC","pPL","coefAge","coefMinT30","coefRain30","pAge","pMinT30","pRain30")
}


#Binomial GLMM for other taxa (that did not converge with a NB model)
ps4<-ps3

#Order
ps4<-droplevels(subset(ps4,ps4$Taxa %in% c("Verrucomicrobiales","Enterobacteriales","Bacillales")))

#Family
ps4<-droplevels(subset(ps4,ps4$Taxa %in% c("Bacteroidaceae","Akkermansiaceae","Enterobacteriaceae","Staphylococcaceae","Marinifilaceae")))

#Genus
ps4<-droplevels(subset(ps4,ps4$Taxa %in% c("Oxalobacter","Alistipes","Butyricimonas","Bacteroides","Akkermansia","[Ruminococcus] gnavus group","Adlercreutzia","Escherichia-Shigella","Staphylococcus","Hungatella")))


#Convert Abundance to binary variable
head(ps4)
unique(ps4$Taxa)
ps4$AbundanceBin<-ps4$Abundance
ps4$AbundanceBin[ps4$AbundanceBin>0]<-1
ps4$AbundanceBin<-as.factor(ps4$AbundanceBin)
table(ps4$AbundanceBin)
 
#Loop
for (i in unique(ps4$Taxa)){
  data<-droplevels(subset(ps4,ps4$Taxa %in% i))
  m1 <- glmer(AbundanceBin ~ ReproState + zAgeYear + zMinT30 + zRain30+ offset(log(NbReads)) + (1|ID) + (1|Unit),data=data,family=binomial)
  shap<-round(as.numeric(testUniformity(simulateResiduals(m1),plot=F)[2]),4)
  var1<-as.data.frame(print(VarCorr(m1),comp=c("Variance","Std.Dev.")))
  varRE<-round(as.numeric(var1[,4]),3)
  sdRE<-round(as.numeric(var1[,5]),3)
  test = glht(m1,linfct=mcp(ReproState="Tukey"))
  pval1<-round(as.numeric(summary(test,test=adjusted(type="none"))$test$pvalues),4)
  coef1<-round(as.numeric(summary(test,test=adjusted(type="none"))$test$coef),4)
  coef<-round(summary(m1)$coef[4:6,1],4)
  drop<-drop1(m1,test="Chisq")
  pval<-round(drop[3:5,4],4)
  df<-c(i,shap,varRE,sdRE,coef1,pval1,coef,pval)
  df2<-rbind(df2,df)
  colnames(df2)<-c("Taxa","pShap","varRE.ID","varRE.Unit","sdRE.ID","sdRE.Unit","coefLC","coefPC","coefPL","pLC","pPC","pPL","coefAge","coefMinT30","coefRain30","pAge","pMinT30","pRain30")
}

#Results
df2<-data.frame(df2)
head(df2)
dim(df2)
length(unique(ps3$Taxa))

#Adjust pvalues for multiple testing
df2$pLC<-as.numeric(as.character(df2$pLC))
df2$pPC<-as.numeric(as.character(df2$pPC))
df2$pPL<-as.numeric(as.character(df2$pPL))
df2$pAge<-as.numeric(as.character(df2$pAge))
df2$pMinT30<-as.numeric(as.character(df2$pMinT30))
df2$pRain30<-as.numeric(as.character(df2$pRain30))
df2$coefLC<-as.numeric(as.character(df2$coefLC))
df2$coefPC<-as.numeric(as.character(df2$coefPC))
df2$coefPL<-as.numeric(as.character(df2$coefPL))
df2$coefAge<-as.numeric(as.character(df2$coefAge))
df2$coefMinT30<-as.numeric(as.character(df2$coefMinT30))
df2$coefRain30<-as.numeric(as.character(df2$coefRain30))

df2$pLC<-p.adjust(df2$pLC,method="fdr")
df2$pPC<-p.adjust(df2$pPC,method="fdr")
df2$pPL<-p.adjust(df2$pPL,method="fdr")
df2$pAge<-p.adjust(df2$pAge,method="fdr")
df2$pMinT30<-p.adjust(df2$pMinT30,method="fdr")
df2$pRain30<-p.adjust(df2$pRain30,method="fdr")

#Adjusted p-values
dim(df2[df2$pLC<=0.05,])
dim(df2[df2$pPC<=0.05,])
dim(df2[df2$pPL<=0.05,])
dim(df2[df2$pMinT30<=0.05,])
dim(df2[df2$pRain30<=0.05,])
dim(df2[df2$pAge<=0.05,])





###################################### 4. FUNCTIONAL ANALYSIS (picrust2) #######################################
#Choose LEVEL 2 or LEVEL 3 dataset
jd<-read.csv("~/Desktop/ko_L2.txt",header=TRUE, sep=",", na.strings="NA")
#jd<-read.csv("~/Desktop/ko_L3.txt",header=TRUE, sep=",", na.strings="NA")
jd[1:5,1:5]
rownames(jd)<-jd$X
jd<-jd[,-1]
jd[1:5,1:5]

## Restrict to female samples
head(meta)
dim(meta) #should be 439 samples
jd<-droplevels(jd[,colnames(jd) %in% meta$SampleID])

#Convert pathway counts to relative abundance
jd[1:5,1:5]
jd1<-data.frame(apply(jd,2,function(x){(x/sum(x))*100}))
jd1[1:5,1:5]#row=pathways,col=samples
dim(jd1)

#Calculate relative abundance of each pathway per sample
jd1bis<-jd1
jd1bis$pathway<-rownames(jd1bis)
jd2<-melt(jd1bis,id.vars="pathway")
head(jd2)
colnames(jd2)<-c("pathway","SampleID","RelativeAbundance")
a1<-sapply(unique(jd2$pathway),function(x) mean(jd2$RelativeAbundance[jd2$pathway==x],na.rm=TRUE))
a2<-data.frame(cbind(Pathway=as.character(unique(jd2$pathway)),RA=a1))
a2$RA<-as.numeric(as.character(a2$RA))
a2[with(a2, order(-RA)), ] 

#Filter to only include pathways >0.1% relative abundance
jd3<-jd1[rownames(jd1) %in% unique(a2$Pathway[a2$RA>=0.1]), ]
jd4<-data.frame(t(jd3))
jd4[1:5,1:5]#rows=samples,columns=pathways
dim(jd4)


#Add metadata
jd4$SampleID<-rownames(jd4)
jd5<-merge(jd4,meta,by="SampleID",all.x=TRUE,all.y=FALSE)
head(jd5)
dim(jd5)


#Z-transformation of numerical variables
stdize=function(x) {(x - mean(x))/sd(x)}
jd5$zAgeYear<-stdize(jd5$AgeYear)
jd5$zRain30<-stdize(jd5$Rain30)
jd5$zMinT30<-stdize(jd5$MinT30)

#Loop
pi<-NULL
for (i in 2:33){ #at LEVEL 2
#for (i in 2:135){ #at LEVEL 3
  m1<-lmer(as.numeric(jd5[,i]) ~ ReproState + zAgeYear + zRain30 + zMinT30 + (1|ID)+ (1|Unit),data=jd5)
  test = glht(m1,linfct=mcp(ReproState="Tukey"))
  coef1<-as.numeric(summary(test,test=adjusted(type="none"))$test$coef)
  pval1<-as.numeric(summary(test,test=adjusted(type="none"))$test$pvalues)
  pval<-drop1(m1,test="Chisq")
  dim(pval)
  pAge<-round(pval[2,4],4)
  pRain30<-round(pval[3,4],4)
  pMinT30<-round(pval[4,4],4)
  coef<-round(coef(summary(m1))[ ,"Estimate"],6)
  res<-c(as.character(i),colnames(jd5)[i],coef1,pval1,coef[4:6],pAge,pRain30,pMinT30)
  pi<-rbind(pi,res)
  colnames(pi)<-c("Col","Taxa","coefLC","coefPC","coefPL","pLC","pPC","pPL","coefAge","coefRain30","coefMinT30","pAge","pRain30","pMinT30")
}
pi<-data.frame(pi)
head(pi)
dim(pi)


#Convert class
pi$pLC<-as.numeric(as.character(pi$pLC))
pi$pPC<-as.numeric(as.character(pi$pPC))
pi$pPL<-as.numeric(as.character(pi$pPL))
pi$pAge<-as.numeric(as.character(pi$pAge))
pi$pRain30<-as.numeric(as.character(pi$pRain30))
pi$pMinT30<-as.numeric(as.character(pi$pMinT30))
pi$coefLC<-as.numeric(as.character(pi$coefLC))
pi$coefPC<-as.numeric(as.character(pi$coefPC))
pi$coefPL<-as.numeric(as.character(pi$coefPL))
pi$coefAge<-as.numeric(as.character(pi$coefAge))
pi$coefRain30<-as.numeric(as.character(pi$coefRain30))
pi$coefMinT30<-as.numeric(as.character(pi$coefMinT30))

#Adjust pvalues for multiple testing
pi$pLC<-p.adjust(pi$pLC,method="fdr")
pi$pPC<-p.adjust(pi$pPC,method="fdr")
pi$pPL<-p.adjust(pi$pPL,method="fdr")
pi$pAge<-p.adjust(pi$pAge,method="fdr")
pi$pRain30<-p.adjust(pi$pRain30,method="fdr")
pi$pMinT30<-p.adjust(pi$pMinT30,method="fdr")

#Adjusted
dim(pi[pi$pLC<=0.05,])
dim(pi[pi$pPC<=0.05,])
dim(pi[pi$pPL<=0.05,])
dim(pi[pi$pAge<=0.05,])
dim(pi[pi$pRain30<=0.05,])
dim(pi[pi$pMinT30<=0.05,])

