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
theme_set (theme_classic(base_size=20))


############## DATASET #################
load("~/Desktop/gelada_SILVA_ASV filter.RData")
gelada_physeq #758 samples - 3295 ASVs
otu_table(gelada_physeq)[1:5,1:5]#rows=ASVs,columns=samples
head(sample_data(gelada_physeq))
summary(sample_data(gelada_physeq))
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
m1<- as.matrix(t(otu_table(gelada_physeq)))
m1[1:5,1:5]#row=samples
taxa<-data.frame(tax_table(gelada_physeq))
taxa$taxa_id<-rownames(taxa)

### clr-transfornmation of the counts ###
#for more information on compositional microbiome analysis, see Gloor et al. 2017 and Gloor's github tutorial (https://github.com/ggloor/CoDa_microbiome_tutorial/wiki/Part-1:-Exploratory-Compositional-PCA-biplot)
pseudocount <- 0.65 #add a pseudo count for zeroes
clr <- t(apply(t(m1)+pseudocount, 2, compositions::clr))
clr[1:5,1:5]#row=samples
dist(clr)#Aitchison distance

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
df1[1:5,c(1:3,760:773)]

#PCA representation
a<-ggplot(df1, aes(x=Rain30,y=PC1)) + geom_point() + ylab(PC1) + xlab("Cumulative rainfall") + theme(axis.title=element_text(size=18,face="bold"))+scale_x_continuous(breaks=seq(0,600,100),limits=c(0,600))
b<-ggplot(df1, aes(x=Rain30,y=PC2)) + geom_point() + ylab(PC1) + xlab("Cumulative rainfall") + theme(axis.title=element_text(size=18,face="bold"))+scale_x_continuous(breaks=seq(0,600,100),limits=c(0,600))
ggarrange(a,b,nrow=1)

a<-ggplot(df1, aes(x=MinT30,y=PC1)) + geom_point() + ylab(PC1) + xlab("Minimum temperature") + theme(axis.title=element_text(size=18,face="bold"))+scale_x_continuous(breaks=seq(5.8,12,2),limits=c(5.8,12))
b<-ggplot(df1, aes(x=MinT30,y=PC2)) + geom_point() + ylab(PC2) + xlab("Minimum temperature") + theme(axis.title=element_text(size=18,face="bold"))+scale_x_continuous(breaks=seq(5.8,12,2),limits=c(5.8,12))
ggarrange(a,b,nrow=1)

a<-ggplot(df1, aes(x=PC1,y=PC2,colour=Sex)) + geom_point()+ xlab(PC1)+ ylab(PC2)+theme(axis.title=element_text(size=18,face="bold"))
b<-ggplot(df1, aes(x=Sex,y=PC1,fill=Sex)) + geom_boxplot()+ ylab(PC1)+xlab("Sex")+theme(axis.title=element_text(size=18,face="bold"))+theme(legend.position="none")
c<-ggplot(df1, aes(x=Sex,y=PC2,fill=Sex)) + geom_boxplot()+ ylab(PC2)+xlab("Sex")+theme(axis.title=element_text(size=18,face="bold"))+theme(legend.position="none")
ggarrange(a,b,c,nrow=1)

a<-ggplot(df1, aes(x=AgeYear,y=PC1)) + geom_point()+ ylab(PC1)+xlab("Age (in years)")+theme(axis.title=element_text(size=18,face="bold")) + scale_x_continuous(breaks=seq(5,30,5),limits=c(4.5,30.2))
b<-ggplot(df1, aes(x=AgeYear,y=PC2)) + geom_point()+ ylab(PC2)+xlab("Age (in years)")+theme(axis.title=element_text(size=18,face="bold")) + scale_x_continuous(breaks=seq(5,30,5),limits=c(4.5,30.2))
ggarrange(a,b,nrow=1)


### Loading scores of taxa on PC1 ###
pca$rotation[1:5,1:5]
t1<-data.frame(pca$rotation)
t1$taxa_id<-rownames(t1)
t1[1:5,1:5]
barplot(sort(t1$PC1,decreasing=TRUE),ylab="Loading PC1")
taxa1<-merge(taxa,t1[c("taxa_id","PC1")],by="taxa_id")
taxa2<-taxa1[with(taxa1,order(-PC1)), ] 
head(taxa2,n=20)


### PERMANOVA test ####
aich<-dist(clr) #Aitchison distance

#Effect of individuals
adonis2(aich ~ ID + NbReads, data = meta, permutations = 10000)

#Effects of predictors
adonis2(aich ~ NbReads + Unit +  Rain30 + MinT30 + Sex + AgeYear, strata=ID, data = meta, permutations = 10000)




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
#Phylum: all phyla converge 

#Class
ps4<-droplevels(subset(ps4,ps4$Taxa %!in% c("Actinobacteria")))

#Order
ps4<-droplevels(subset(ps4,ps4$Taxa %!in% c("Enterobacteriales","Bifidobacteriales","Bacillales")))

#Family
ps4<-droplevels(subset(ps4,ps4$Taxa %!in% c("Bacteroidaceae","Enterobacteriaceae","Bifidobacteriaceae","Staphylococcaceae")))

#Genus
ps4<-droplevels(subset(ps4,ps4$Taxa %!in% c("Bacteroides","Alistipes","Staphylococcus","[Ruminococcus] gnavus group","Bifidobacterium","Escherichia-Shigella")))


### LOOP ###
#Fit a NB GLMM for each taxa
df2<-NULL
for (i in unique(ps4$Taxa)){
  data<-droplevels(subset(ps4,ps4$Taxa %in% i))
  m1 <- glmer.nb (Abundance ~ Sex + zAgeYear + zMinT30 + zRain30 + offset(log(NbReads)) + (1|ID) + (1|Unit),data=data)
  shap<-round(as.numeric(testUniformity(simulateResiduals(m1),plot=F)[2]),4)
  var1<-as.data.frame(print(VarCorr(m1),comp=c("Variance","Std.Dev.")))
  varRE<-round(as.numeric(var1[,4]),3)
  sdRE<-round(as.numeric(var1[,5]),3)
  coef<-round(summary(m1)$coef[2:5,1],4)
  drop<-drop1(m1,test="Chisq")
  pval<-round(drop[2:5,4],4)
  df<-c(i,shap,varRE,sdRE,coef,pval)
  df2<-rbind(df2,df)
  colnames(df2)<-c("Taxa","pShap","varRE.ID","varRE.Unit","sdRE.ID","sdRE.Unit","coefSexM","coefAge","coefMinT30","coefRain30","pSex","pAge","pMinT30","pRain30")
}


#Binomial GLMM for other taxa (that did not converge with a NB model)
ps4<-ps3

#Class
ps4<-droplevels(subset(ps4,ps4$Taxa %in% c("Actinobacteria")))

#Order
ps4<-droplevels(subset(ps4,ps4$Taxa %in% c("Enterobacteriales","Bifidobacteriales","Bacillales")))

#Family
ps4<-droplevels(subset(ps4,ps4$Taxa %in% c("Bacteroidaceae","Enterobacteriaceae","Bifidobacteriaceae","Staphylococcaceae")))

#Genus
ps4<-droplevels(subset(ps4,ps4$Taxa %in% c("Bacteroides","Alistipes","Staphylococcus","[Ruminococcus] gnavus group","Bifidobacterium","Escherichia-Shigella")))

#Convert Abundance to binary variable
ps4$AbundanceBin<-ps4$Abundance
ps4$AbundanceBin[ps4$AbundanceBin>0]<-1
ps4$AbundanceBin<-as.factor(ps4$AbundanceBin)

#Loop
for (i in unique(ps4$Taxa)){
  data<-droplevels(subset(ps4,ps4$Taxa %in% i))
  m1 <- glmer(AbundanceBin ~ Sex + zAgeYear + zMinT30 + zRain30 + offset(log(NbReads)) + (1|ID) + (1|Unit),data=data,family=binomial)
  shap<-round(as.numeric(testUniformity(simulateResiduals(m1),plot=F)[2]),4)
  var1<-as.data.frame(print(VarCorr(m1),comp=c("Variance","Std.Dev.")))
  varRE<-round(as.numeric(var1[,4]),3)
  sdRE<-round(as.numeric(var1[,5]),3)
  coef<-round(summary(m1)$coef[2:5,1],4)
  drop<-drop1(m1,test="Chisq")
  pval<-round(drop[2:5,4],4)
  df<-c(i,shap,varRE,sdRE,coef,pval)
  df2<-rbind(df2,df)
  colnames(df2)<-c("Taxa","pShap","varRE.ID","varRE.Unit","sdRE.ID","sdRE.Unit","coefSexM","coefAge","coefMinT30","coefRain30","pSex","pAge","pMinT30","pRain30")
}


#Results
df2<-data.frame(df2)
head(df2)
dim(df2)
length(unique(ps3$Taxa))

#Adjust pvalues for multiple testing
df2$pSex<-as.numeric(as.character(df2$pSex))
df2$pAge<-as.numeric(as.character(df2$pAge))
df2$pMinT30<-as.numeric(as.character(df2$pMinT30))
df2$pRain30<-as.numeric(as.character(df2$pRain30))
df2$coefSex<-as.numeric(as.character(df2$coefSex))
df2$coefAge<-as.numeric(as.character(df2$coefAge))
df2$coefMinT30<-as.numeric(as.character(df2$coefMinT30))
df2$coefRain30<-as.numeric(as.character(df2$coefRain30))

df2$pSex<-p.adjust(df2$pSex,method="fdr")
df2$pAge<-p.adjust(df2$pAge,method="fdr")
df2$pMinT30<-p.adjust(df2$pMinT30,method="fdr")
df2$pRain30<-p.adjust(df2$pRain30,method="fdr")

#Adjusted p-values
dim(df2[df2$pSex<=0.05,])
dim(df2[df2$pAge<=0.05,])
dim(df2[df2$pMinT30<=0.05,])
dim(df2[df2$pRain30<=0.05,])


### FIGURES ###

#Significant estimates of rainfall
head(df2)
rain<-df2[,c("Taxa","pRain30","coefRain30")]
rain<-droplevels(rain[rain$pRain30<=0.05,])
rain<-rain[with(rain,order(coefRain30)), ] 
rain
rain$col<-ifelse(rain$coefRain30<0,"increase in dry periods","increase in wet periods")
rain$Taxa<- factor(rain$Taxa, levels = rain$Taxa)

ggplot(rain, aes(x = Taxa, y = coefRain30,fill=col))+ xlab("taxa") + ylab("Estimate of cumulative rainfall") +
  geom_bar(stat = 'identity', position = "identity") + coord_flip() +
  scale_fill_manual("",breaks=c("increase in wet periods","increase in dry periods"),values = c("increase in dry periods" = "darkgoldenrod1","increase in wet periods" = "cornflowerblue"))+
  theme(legend.position = "bottom",legend.background=element_blank(),legend.text = element_text( size = 10, face = "bold"),axis.text=element_text(size=13,face="bold"),axis.title=element_text(size=16,face="bold"))

#Significant estimates of temperature
temp<-df2[,c("Taxa","pMinT30","coefMinT30")]
temp<-droplevels(temp[temp$pMinT30<=0.05,])
temp<-temp[with(temp,order(coefMinT30)), ] 
temp
temp$col<-ifelse(temp$coefMinT30<0,"increase in cold periods","increase in hot periods")
temp$Taxa<- factor(temp$Taxa, levels = temp$Taxa)

ggplot(temp, aes(x = Taxa, y = coefMinT30,fill=col))+ xlab("Taxa")+ ylab("Estimate of minimum temperature")+
  geom_bar(stat = 'identity', position = "identity") + coord_flip() +
  scale_fill_manual("",breaks=c("increase in cold periods","increase in hot periods"),values = c("increase in hot periods" = "darkgoldenrod1","increase in cold periods" = "cornflowerblue"))+
  theme(legend.position = "bottom",legend.background=element_blank(),legend.text = element_text( size = 10, face = "bold"),axis.text=element_text(size=13,face="bold"),axis.title=element_text(size=16,face="bold"))




#RAW COUNT DATA WITH RAINFALL
fam = tax_glom(gelada_physeq,"Family",NArm=FALSE)
fam1<-psmelt(fam)
fam1$RelativeAbundance<-(as.numeric(as.character(fam1$Abundance))*100)/as.numeric(as.character(fam1$NbReads))

gen = tax_glom(gelada_physeq,"Genus",NArm=FALSE)
gen1<-psmelt(gen)
gen1$RelativeAbundance<-(as.numeric(as.character(gen1$Abundance))*100)/as.numeric(as.character(gen1$NbReads))

df1<-droplevels(subset(fam1,fam1$Family %in% "Prevotellaceae"))
a1<-ggplot(df1, aes(x=Rain30,y=RelativeAbundance+0.01)) + geom_point() + xlim(0,600) + ggtitle("Prevotellaceae")+ xlab("Cumulative rainfall") + ylab("Relative abundance (%)")+ geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=18,face="bold"))

df2<-droplevels(subset(fam1,fam1$Family %in% "Lachnospiraceae"))
a2<-ggplot(df2, aes(x=Rain30,y=RelativeAbundance+0.01)) + geom_point() + xlim(0,600) + ggtitle("Lachnospiraceae")+ xlab("Cumulative rainfall")+ ylab("Relative abundance (%)") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=18,face="bold"))

df3<-droplevels(subset(fam1,fam1$Family %in% "Fibrobacteraceae"))
a3<-ggplot(df3, aes(x=Rain30,y=RelativeAbundance+0.01)) + geom_point() + xlim(0,600) + ggtitle("Fibrobacteraceae")+ xlab("Cumulative rainfall")+ ylab("Relative abundance (%)") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=18,face="bold"))

df4<-droplevels(subset(gen1,gen1$Genus %in% "Succinivibrio"))
a4<-ggplot(df4, aes(x=Rain30,y=RelativeAbundance+0.01)) + geom_point() + xlim(0,600) + ggtitle("Succinivibrio")+ xlab("Cumulative rainfall")+ ylab("Relative abundance (%)") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=18,face="bold"))

df5<-droplevels(subset(gen1,gen1$Genus %in% "Methanobrevibacter"))
a5<-ggplot(df5, aes(x=Rain30,y=RelativeAbundance+0.01)) + geom_point() + xlim(0,600) + ggtitle("Methanobrevibacter")+ xlab("Cumulative rainfall")+ ylab("Relative abundance (%)") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=18,face="bold"))

df6<-droplevels(subset(fam1,fam1$Family %in% "RFP12"))
a6<-ggplot(df6, aes(x=Rain30,y=RelativeAbundance+0.01)) + geom_point() + ggtitle("RFP12")+ xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + geom_smooth(method='lm')+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=18,face="bold"))

ggarrange(a1,a2,a3,a4,a5,a6,labels=c("A","","","B","",""),ncol=3,nrow=2)




###################################### 4. FUNCTIONAL ANALYSIS (picrust2) #######################################
#Choose LEVEL 2 or LEVEL 3 dataset
jd<-read.csv("~/Desktop/ko_L2.txt",header=TRUE, sep=",", na.strings="NA")
#jd<-read.csv("~/Desktop/ko_L3.txt",header=TRUE, sep=",", na.strings="NA")
jd[1:5,1:5]
rownames(jd)<-jd$X
jd<-jd[,-1]
jd[1:5,1:5]


#Convert pathway counts to relative abundance
jd[1:5,1:5]
jd1<-data.frame(apply(jd,2,function(x){(x/sum(x))*100}))
jd1[1:5,1:5]#row=pathways,col=samples

#Calculate relative abundance of each pathway per sample
jd1bis<-jd1
jd1bis$pathway<-rownames(jd1bis)
jd2<-melt(jd1bis,id.vars="pathway")
colnames(jd2)<-c("pathway","SampleID","RelativeAbundance")
a1<-sapply(unique(jd2$pathway),function(x) mean(jd2$RelativeAbundance[jd2$pathway==x],na.rm=TRUE))
a2<-data.frame(cbind(Pathway=as.character(unique(jd2$pathway)),RA=a1))
a2$RA<-as.numeric(as.character(a2$RA))
a2[with(a2, order(-RA)), ] 
dim(a2)

#Filter to only include pathways >0.1% relative abundance
jd3<-jd1[rownames(jd1) %in% unique(a2$Pathway[a2$RA>=0.1]), ]
jd4<-data.frame(t(jd3))
jd4[1:5,1:5]#rows=samples,columns=pathways
dim(jd4)

#Add metadata
jd4$SampleID<-rownames(jd4)
jd5<-merge(jd4,meta,by="SampleID",all.x=TRUE,all.y=FALSE)
jd5<-droplevels(subset(jd5,jd5$SampleID %!in% NA))
head(jd5)
dim(jd5)

#Z-transformation of numerical variables
stdize=function(x) {(x - mean(x))/sd(x)}
jd5$zAgeYear<-stdize(jd5$AgeYear)
jd5$zRain30<-stdize(jd5$Rain30)
jd5$zMinT30<-stdize(jd5$MinT30)


# LMM #
names(jd5)
pi<-NULL
for (i in 2:33){ #for LEVEL 2
#for (i in 2:135){ #for LEVEL 3
  m1<-lmer(as.numeric(jd5[,i]) ~ Sex + zAgeYear+ zRain30 + zMinT30 + (1|ID)+ (1|Unit),data=jd5)
  shap<-round(as.numeric(shapiro.test(resid(m1))[2]),4)
  pval<-drop1(m1,test="Chisq")
  pSex<-round(pval[2,4],4)
  pAge<-round(pval[3,4],4)
  pRain30<-round(pval[4,4],4)
  pMinT30<-round(pval[5,4],4)
  coef1<-round(coef(summary(m1))[ ,"Estimate"],6)
  res<-c(as.character(i),colnames(jd5)[i],shap,pSex,pAge,pRain30,pMinT30,coef1[2:5])
  pi<-rbind(pi,res)
  colnames(pi)<-c("Col","Pathway","pShap","pSex","pAge","pRain30","pMinT30","coefSex","coefAge","coefRain30","coefMinT30")
}
pi<-data.frame(pi)
head(pi)
dim(pi)


#Convert class
pi$pShap<-as.numeric(as.character(pi$pShap))
pi$pSex<-as.numeric(as.character(pi$pSex))
pi$pAge<-as.numeric(as.character(pi$pAge))
pi$pRain30<-as.numeric(as.character(pi$pRain30))
pi$pMinT30<-as.numeric(as.character(pi$pMinT30))
pi$coefSex<-as.numeric(as.character(pi$coefSex))
pi$coefAge<-as.numeric(as.character(pi$coefAge))
pi$coefRain30<-as.numeric(as.character(pi$coefRain30))
pi$coefMinT30<-as.numeric(as.character(pi$coefMinT30))

#Adjust pvalues for multiple testing
pi$pSex<-p.adjust(pi$pSex,method="fdr")
pi$pAge<-p.adjust(pi$pAge,method="fdr")
pi$pRain30<-p.adjust(pi$pRain30,method="fdr")
pi$pMinT30<-p.adjust(pi$pMinT30,method="fdr")

#Signiifcant pathways with predictors
dim(pi[pi$pRain30<=0.05,])
dim(pi[pi$pMinT30<=0.05,])
dim(pi[pi$pSex<=0.05,])
dim(pi[pi$pAge<=0.05,])



### FIGURES ###
pi$Pathway<-gsub('\\.', ' ',pi$Pathway)

#Significant estimates of rainfall
rain<-pi[,c("Pathway","pRain30","coefRain30")]
rain<-droplevels(rain[rain$pRain30<=0.05,])
rain<-rain[with(rain,order(coefRain30)), ] 
rain

rain$col<-ifelse(rain$coefRain30<0,"increase in dry periods","increase in wet periods")
rain$Pathway<- factor(rain$Pathway, levels = rain$Pathway)

ggplot(rain, aes(x = Pathway, y = coefRain30,fill=col))+ xlab("Pathway")+ ylab("Estimate of cumulative rainfall")+
  geom_bar(stat = 'identity', position = "identity") + coord_flip()+
  scale_fill_manual("",breaks=c("increase in wet periods","increase in dry periods"),values = c("increase in dry periods" = "darkgoldenrod1","increase in wet periods" = "cornflowerblue"))+ 
  theme(legend.position = "bottom",legend.background=element_blank(),legend.text = element_text( size = 10, face = "bold"),axis.text=element_text(size=13,face="bold"),axis.title=element_text(size=16,face="bold"))



#Significant estimates of temperature
temp<-pi[,c("Pathway","pMinT30","coefMinT30")]
temp<-droplevels(temp[temp$pMinT30<=0.05,])
temp<-temp[with(temp,order(coefMinT30)), ] 
temp
temp$col<-ifelse(temp$coefMinT30<0,"increase in cold periods","increase in hot periods")
temp$Pathway<- factor(temp$Pathway, levels = temp$Pathway)

ggplot(temp, aes(x = Pathway, y = coefMinT30,fill=col))+ xlab("Pathway")+ ylab("Estimate of minimum temperature")+
  geom_bar(stat = 'identity', position = "identity") + coord_flip() +
  scale_fill_manual("",breaks=c("increase in cold periods","increase in hot periods"),values = c("increase in hot periods" = "darkgoldenrod1","increase in cold periods" = "cornflowerblue"))+
  theme(legend.position = "bottom",legend.background=element_blank(),legend.text = element_text( size = 10, face = "bold"),axis.text=element_text(size=13,face="bold"),axis.title=element_text(size=16,face="bold"))


#RAW COUNT DATA WITH RAINFALL
#LEVEL 2
jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Energy.Metabolism")]
a1<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Energy Metabolism") + geom_smooth(method='lm')+scale_y_log10()+ theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=14,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Amino.Acid.Metabolism")]
a2<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Amino Acid Metabolism") + geom_smooth(method='lm')+scale_y_log10()+ theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=14,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Lipid.Metabolism")]
a3<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Lipid Metabolism") + geom_smooth(method='lm')+scale_y_log10()+  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=14,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Membrane.Transport")]
a4<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Membrane Transport") + geom_smooth(method='lm')+scale_y_log10()+ theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=14,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Replication.and.Repair")]
a5<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Replication and Repair") + geom_smooth(method='lm')+scale_y_log10()+ theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=14,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Cell.Motility")]
a6<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Cell Motility") + geom_smooth(method='lm')+scale_y_log10()+ theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=14,face="bold"))

ggarrange(a1,a2,a3,a4,a5,a6,labels=c("A","","","B","",""),ncol=3,nrow=2)



#LEVEL 3
jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Oxidative.phosphorylation")]
a1<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Oxidative phosphorylation") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Citrate.cycle..TCA.cycle.")]
a2<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Citrate cycle") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Fatty.acid.biosynthesis")]
a3<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Fatty acid biosynthesis") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Carbon.fixation.pathways.in.prokaryotes")]
a4<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Carbon fixation pathways 
in prokaryotes") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Methane.metabolism")]
a5<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Methane metabolism") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Transporters")]
a6<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Membrane Transporters") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Starch.and.sucrose.metabolism")]
a7<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall")+ ylab("Relative abundance (%)e") + ggtitle("Starch and sucrose
metabolism") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"))

jd5$Pathway<-jd5[,which(colnames(jd5) %in% "Fructose.and.mannose.metabolism")]
a8<-ggplot(jd5, aes(x=Rain30,y=Pathway)) + geom_point()  + xlab("Cumulative rainfall") + ylab("Relative abundance (%)") + ggtitle("Fructose and mannose 
metabolism") + geom_smooth(method='lm')+scale_y_log10()+theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"),plot.title = element_text(size=16,face="bold"))

ggarrange(a1,a2,a3,a4,a5,a6,a7,a8,labels=c("A","","","","","B","","","",""),ncol=4,nrow=2)




