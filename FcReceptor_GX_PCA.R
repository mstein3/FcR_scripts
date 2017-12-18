### Plotting Gene expression PCA data for Fc receptor manuscript

# Use the new processed data: 
nt.pheno<-read.table("Pheno_Rready_nullCD_7M.txt", sep="\t", header=TRUE, fill=TRUE)
nt.gx<-read.table("VoomTMM_NullCD_Elist_RMPool.txt", sep="\t")
# Move to a new directory: 
setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_Cyto2/Fcreceptor2/May17_GX_PCA/")

# Add fc.pca.color and treat.shape 
dgreen<-c("#1B9E77") # GG, rgb: 27, 158, 119
dorange<-c("#D95F02") # AA,  rgb: 271, 95, 2
dpurple<-c("#7570B3") #AG,  rgb: 117, 112, 179 

nt.pheno2$fc.pca.color<-ifelse(nt.pheno2$FcConsensus_letter=="AA", dorange, ifelse(nt.pheno2$FcConsensus_letter=="AG", dpurple, dgreen))
nt.pheno2$treat.shape<-ifelse(nt.pheno2$treat=="null", 2, ifelse(nt.pheno2$treat=="CD", 19, 17))

#### 0.5 Also remove one outlier
nt.pheno3<-nt.pheno2[!outlier.to.rm,]
nt.gx3<-nt.gx2[,!outlier.to.rm]

write.table(nt.gx3,"ImHTGX_CDNull_nomissingFc.txt", sep="\t", quote=FALSE)
write.table(nt.pheno3, "Pheno_ImHTGX_CDNull_nomissingFc.txt", sep="\t", quote=FALSE)

###### 1. PCA plots with null + CD3, all 3 genotypes: ==========
## Run PCA:
pheno.file<-nt.pheno3
v.gx<-nt.gx3

# Set up covariates: 
treatment<-as.factor(pheno.file$treat)
fc.geno<-as.numeric(pheno.file$FcConsensus_num)
fc.geno.letter<-as.factor(pheno.file$FcConsensus_letter)
sex<-as.factor(pheno.file$sex)
age<-as.numeric(pheno.file$Age)
extbatch<-as.factor(pheno.file$Extract_Batch)
atopy<-as.factor(pheno.file$Allergy_stat_2isyes)
asthma<-as.factor(pheno.file$Asthma_stat)
RIN<-as.numeric(pheno.file$RIN_final)
adapter_num<-as.factor(pheno.file$adapter_num)
qpcr_conc<-as.numeric(pheno.file$qPCR_ng_ul)
bioA_conc<-as.numeric(pheno.file$BioA_ng_ul)
library_tech<-as.factor(pheno.file$library_tech)
num_reads<-as.numeric(pheno.file$Total_Map)
all.pool.code<-as.factor(pheno.file$Comb_pool_code)
orig.pool.code<-as.factor(pheno.file$Pool_1.0)
all.lane.code<-as.factor(pheno.file$Comb_lane_code_allruns)
all.fc.code<-as.factor(pheno.file$Comb_flowcell_code_allruns)
indiv<-as.factor(pheno.file$FINDIV)

pc_list=c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
gx.covars<-cbind.data.frame(treatment, fc.geno, fc.geno.letter, sex, age, extbatch, atopy, asthma, RIN, adapter_num, qpcr_conc, bioA_conc, library_tech, num_reads, all.pool.code, orig.pool.code, all.lane.code, all.fc.code)

sum.PC <- prcomp(t(v.gx), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:20] 

# To see if covariates are correlated with a PC (looking at PC1-50)
pval.pca1=matrix(ncol=ncol(gx.covars), nrow=20)
rownames(pval.pca1)=pc_list
colnames(pval.pca1)=colnames(gx.covars)

for(j in 1:ncol(gx.covars))
{
  for(i in 1:length(pc_list))
  {
    data1= lm(sum.PC$x[,i]~gx.covars[,j])
    pval.pca1[i,j]=anova(data1)$'Pr(>F)'[1]
  }
}

setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/CDNull_Processing")
write.table(pval.pca1, "20PCs_vs_covar_FcReceptor_TMMrmPool_rm1outlier.txt", sep="\t", quote=FALSE)

# Round sumsum: 
# Modify % variance to show 1 significant digit: 
sumsum.round<-NULL
for(i in 1:length(sumsum$importance[2,])){
  sumsum.round[i]<-round((sumsum$importance[2,i]*100), digits=1)
}

setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_Cyto2/Fcreceptor2/May17_GX_PCA/")
# Plot files: 
pdf("ImHTGX_CDnull_FCrecpetorPCA.pdf")
plot(sum.PC$x[,1], sum.PC$x[,2],pch=pheno.file$treat.shape, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 

plot(sum.PC$x[,1], sum.PC$x[,2],pch=pheno.file$treat.shape, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 
text(sum.PC$x[,1],sum.PC$x[,2], fc.geno.letter, cex = .5, pos=3) 
dev.off()


###### 2. PCA plots with CD3 only, all 3 genotypes: ===============
# And outlier removed: 
cd.samp<-nt.pheno3$treat=="CD"
pheno.119<-nt.pheno3[cd.samp,]
v.119<-nt.gx3[,cd.samp]

pheno.file<-pheno.119
v.gx<-v.119

# Set up covariates: 
fc.geno<-as.numeric(pheno.file$FcConsensus_num)
fc.geno.letter<-as.factor(pheno.file$FcConsensus_letter)
sex<-as.factor(pheno.file$sex)
age<-as.numeric(pheno.file$Age)
extbatch<-as.factor(pheno.file$Extract_Batch)
atopy<-as.factor(pheno.file$Allergy_stat_2isyes)
asthma<-as.factor(pheno.file$Asthma_stat)
RIN<-as.numeric(pheno.file$RIN_final)
adapter_num<-as.factor(pheno.file$adapter_num)
qpcr_conc<-as.numeric(pheno.file$qPCR_ng_ul)
bioA_conc<-as.numeric(pheno.file$BioA_ng_ul)
library_tech<-as.factor(pheno.file$library_tech)
num_reads<-as.numeric(pheno.file$Total_Map)
all.pool.code<-as.factor(pheno.file$Comb_pool_code)
orig.pool.code<-as.factor(pheno.file$Pool_1.0)
all.lane.code<-as.factor(pheno.file$Comb_lane_code_allruns)
all.fc.code<-as.factor(pheno.file$Comb_flowcell_code_allruns)

pc_list=c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
gx.covars<-cbind.data.frame(fc.geno, fc.geno.letter, sex, age, extbatch, atopy, asthma, RIN, adapter_num, qpcr_conc, bioA_conc, library_tech, num_reads, all.pool.code, orig.pool.code, all.lane.code, all.fc.code)

sum.PC <- prcomp(t(v.gx), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:20] 

# To see if covariates are correlated with a PC (looking at PC1-20)
pval.pca1=matrix(ncol=ncol(gx.covars), nrow=20)
rownames(pval.pca1)=pc_list
colnames(pval.pca1)=colnames(gx.covars)

for(j in 1:ncol(gx.covars))
{
  for(i in 1:length(pc_list))
  {
    data1= lm(sum.PC$x[,i]~gx.covars[,j])
    pval.pca1[i,j]=anova(data1)$'Pr(>F)'[1]
  }
}

setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/CDNull_Processing")
write.table(pval.pca1, "20PCs_vs_covar_FcReceptor_CDonly_TMMrmPool.txt", sep="\t", quote=FALSE)

# Round sumsum: 
# Modify % variance to show 1 significant digit: 
sumsum.round<-NULL
for(i in 1:length(sumsum$importance[2,])){
  sumsum.round[i]<-round((sumsum$importance[2,i]*100), digits=1)
}

setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_Cyto2/Fcreceptor2/May17_GX_PCA/")
# Plot files: 
pdf("ImHTGX_PCA_CDsamp_TMMnoPool_FCrecpetor.pdf")
plot(sum.PC$x[,1], sum.PC$x[,2],pch=19, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 

plot(sum.PC$x[,1], sum.PC$x[,2],pch=19, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 
text(sum.PC$x[,1],sum.PC$x[,2], indiv, cex = .5, pos=2) 
#dev.off()

# Make this a ggplot2 as well: 
library(ggplot2)
# Make a file: 
boxplot.gg<-cbind.data.frame(sum.PC$x[,1], fc.geno.letter, pheno.file$fc.pca.color)
colnames.gg<-c("PC1", "fc.genotype", "fc.pca.color")
colnames(boxplot.gg)<-colnames.gg

pc.boxplot<-ggplot(boxplot.gg, aes(x=boxplot.gg$fc.genotype, y=boxplot.gg$PC1, fill=boxplot.gg$fc.pca.color)) + geom_boxplot()+
  guides(fill=FALSE)+ ggtitle("Colored by by FCGR2A genotype \n ") +
  theme(axis.title.x = element_blank(), axis.text=element_text(size=16, face="bold", colour="black"), axis.title=element_text(size=16, face="bold"), plot.title= element_text(lineheight = 0.8, face="bold", size=16))+ 
  scale_x_discrete(limits=c("GG","AG","AA"), labels=c("GG", "AG", "AA")) +
  scale_fill_manual(values=c(dgreen, dpurple, dorange))+
  geom_point(position= position_jitter(width=0.03), size=2)+ ylab(paste("PC 1 -",(sumsum.round[1]),"% of variance", sep=" "))
pc.boxplot
dev.off()



###### 3. PCA plots with CD3 only, AG and GG genotypes only: ===============
# Remove AA's from the anti-CD3 gene expression data. 
aa.genos<-pheno.119$FcConsensus_letter=="AA"
pheno.103<-pheno.119[!aa.genos, ]
v.103<-v.119[,!aa.genos]

pheno.file<-pheno.103
v.gx<-v.103

# Set up covariates: 
fc.geno<-as.numeric(pheno.file$FcConsensus_num)
fc.geno.letter<-as.factor(pheno.file$FcConsensus_letter)
sex<-as.factor(pheno.file$sex)
age<-as.numeric(pheno.file$Age)
extbatch<-as.factor(pheno.file$Extract_Batch)
atopy<-as.factor(pheno.file$Allergy_stat_2isyes)
asthma<-as.factor(pheno.file$Asthma_stat)
RIN<-as.numeric(pheno.file$RIN_final)
adapter_num<-as.factor(pheno.file$adapter_num)
qpcr_conc<-as.numeric(pheno.file$qPCR_ng_ul)
bioA_conc<-as.numeric(pheno.file$BioA_ng_ul)
library_tech<-as.factor(pheno.file$library_tech)
num_reads<-as.numeric(pheno.file$Total_Map)
all.pool.code<-as.factor(pheno.file$Comb_pool_code)
orig.pool.code<-as.factor(pheno.file$Pool_1.0)
all.lane.code<-as.factor(pheno.file$Comb_lane_code_allruns)
all.fc.code<-as.factor(pheno.file$Comb_flowcell_code_allruns)

pc_list=c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")
gx.covars<-cbind.data.frame(fc.geno, fc.geno.letter, sex, age, extbatch, atopy, asthma, RIN, adapter_num, qpcr_conc, bioA_conc, library_tech, num_reads, all.pool.code, orig.pool.code, all.lane.code, all.fc.code)

sum.PC <- prcomp(t(v.gx), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:20] 

# Add information to PCA_initial excel file in R_processing 

# To see if covariates are correlated with a PC (looking at PC1-20)
pval.pca1=matrix(ncol=ncol(gx.covars), nrow=20)
rownames(pval.pca1)=pc_list
colnames(pval.pca1)=colnames(gx.covars)

for(j in 1:ncol(gx.covars))
{
  for(i in 1:length(pc_list))
  {
    data1= lm(sum.PC$x[,i]~gx.covars[,j])
    pval.pca1[i,j]=anova(data1)$'Pr(>F)'[1]
  }
}

setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/CDNull_Processing")
write.table(pval.pca1, "20PCs_vs_covar_FcReceptor_CDonly_AGvsGG.txt", sep="\t", quote=FALSE)

# Round sumsum: 
# Modify % variance to show 1 significant digit: 
sumsum.round<-NULL
for(i in 1:length(sumsum$importance[2,])){
  sumsum.round[i]<-round((sumsum$importance[2,i]*100), digits=1)
}

setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_Cyto2/Fcreceptor2/May17_GX_PCA/")

# Plot PCA plots 
pdf("ImHTGX_PCA_CDsamp_TMMnoPool_FCrecpetor_noAAs.pdf")
plot(sum.PC$x[,1], sum.PC$x[,2],pch=19, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 
plot(sum.PC$x[,1], sum.PC$x[,2],pch=19, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 
text(sum.PC$x[,1],sum.PC$x[,2], indiv, cex = .5, pos=2) 


# Make this a ggplot2 as well: 
library(ggplot2)
# Make a file: 
boxplot.gg<-cbind.data.frame(sum.PC$x[,1], sum.PC$x[,2], factor(pheno.file$FcConsensus_letter), factor(pheno.file$fc.pca.color))
colnames.gg<-c("PC1","PC2", "fc.genotype", "fc.pca.color")
colnames(boxplot.gg)<-colnames.gg

pc.boxplot<-ggplot(boxplot.gg, aes(x=boxplot.gg$fc.genotype, y=boxplot.gg$PC1, fill=boxplot.gg$fc.pca.color)) + geom_boxplot()+
  guides(fill=FALSE)+ ggtitle("Colored by by FCGR2A genotype \n ") +
  theme(axis.title.x = element_blank(), axis.text=element_text(size=16, face="bold", colour="black"), axis.title=element_text(size=16, face="bold"), plot.title= element_text(lineheight = 0.8, face="bold", size=16))+ 
  scale_x_discrete(limits=c("GG","AG"), labels=c("GG", "AG")) +
  scale_fill_manual(values=c(dgreen, dpurple))+
  geom_point(position= position_jitter(width=0.03), size=2)+ ylab(paste("PC 1 -",(sumsum.round[1]),"% of variance", sep=" "))

pc2.boxplot<-ggplot(boxplot.gg, aes(x=boxplot.gg$fc.genotype, y=boxplot.gg$PC2, fill=boxplot.gg$fc.pca.color)) + geom_boxplot()+
  guides(fill=FALSE)+ ggtitle("Colored by by FCGR2A genotype \n ") +
  theme(axis.title.x = element_blank(), axis.text=element_text(size=16, face="bold", colour="black"), axis.title=element_text(size=16, face="bold"), plot.title= element_text(lineheight = 0.8, face="bold", size=16))+ 
  scale_x_discrete(limits=c("GG","AG"), labels=c("GG", "AG")) +
  scale_fill_manual(values=c(dgreen,dpurple))+
  geom_point(position= position_jitter(width=0.03), size=2)+ ylab(paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "))

pc.boxplot
pc2.boxplot
dev.off()