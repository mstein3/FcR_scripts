### Processing untreated and CD3 RNAseq data
# Read in gene count files and phenotype information- subset by treatment, samples with 7M reads or more.  
setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/PreQC_GeneCount_pheno_files/")

nt.pheno<-read.table("Pheno_Rready_nullCD_7M.txt", sep="\t", header=TRUE, fill=TRUE)
nt.gc<-read.table("RawGeneCounts_NullCD_RunsCombined_GeneNames_noQC.txt", sep="\t")

## Add colors: 
dgreen<-c("#1B9E77") # GG
#rgb: 27, 158, 119
dorange<-c("#D95F02") # AA
# rgb: 271, 95, 2
dpurple<-c("#7570B3") #AG
# rgb: 117, 112, 179 
nt.pheno$fc.pca.color<-ifelse(nt.pheno$Fc_letter=="AA", dorange, ifelse(nt.pheno$Fc_letter=="AG", dpurple, dgreen))
nt.pheno$treat.shape<-ifelse(nt.pheno$treat=="null", 2, ifelse(nt.pheno$treat=="CD", 19, 17))

######### PART A: Processing ===========================
##### 1. Remove genes with <1 CPM in XXX samples: ===========
# Convert to CPM: 
counts1 = nt.gc+1
CPM.nt.gc = apply(counts1, 2, function(x) (x*1000000)/sum(x))

less.one<-NULL
summary.low.count<-matrix(nrow=nrow(CPM.nt.gc), ncol=1)
for(i in 1:nrow(CPM.nt.gc)){
  less.one<-NULL
  less.one<-which(CPM.nt.gc[i,]<1)
  summary.low.count[i,]<-length(less.one)
  rm(less.one)
}
rownames(summary.low.count)<-rownames(nt.gc)

# Determine which samples have count <1 in at least 20 samples: 
keep.genes<-summary.low.count[,1]<20
table(keep.genes)
#FALSE  TRUE 
# 8199 10328 

# Keep only the genes detected as expressed (removed genes with <1 CPM in at 20+ samples) 
nt.gc2<-nt.gc[keep.genes,]

##### 2. Convert from ENSG to gene names : ===========
ensg.auto2 <-read.table("ENSG_key_ProteinCodingAutosomes_18K_overlap.txt", sep="\t")

ensg.nt.intersect<-intersect(rownames(nt.gc2), ensg.auto2[,4])
ensg.nt.match<-match(ensg.nt.intersect, ensg.auto2[,4])
ensg.auto.nt<-ensg.auto2[ensg.nt.match,]
length(unique(ensg.auto.nt$V6))
# [1] 10328 
# Add gene names: 
nt.gc3<-nt.gc2
rownames(nt.gc3)<-ensg.auto.nt[,6]

setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/CDNull_Processing/")
write.table(nt.gc3, "RawGC_DetectGenes_NullCD_GeneNames.txt", sep="\t", quote=FALSE)

######### PART B: Processing for DE ===========================
###### 5. Normalize =============================
## Normalization of counts: 
library(limma)
library(edgeR)

# To calculate the TMM normalization factors, create a DGElist using the edgeR package: 
dge.nt <- DGEList(counts=nt.gc3)
dge.nt <- calcNormFactors(dge.nt)

# Apply the voom transformation 
v.nt <- voom(dge.nt ,design=NULL,plot=TRUE)
saveRDS(v.nt, "VoomTMM_NullCD_Elist_NoAdj_52417.rds")

######### 6. PCA this data ==============
pheno.file<-nt.pheno
v.gx<-v.nt$E

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

pc_list=c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30", "PC31", "PC32", "PC33", "PC34", "PC35", "PC36", "PC37", "PC38", "PC39", "PC40", "PC41", "PC42", "PC43", "PC44", "PC45", "PC46", "PC47", "PC48", "PC49", "PC50")
gx.covars<-cbind.data.frame(treatment, fc.geno, fc.geno.letter, sex, age, extbatch, atopy, asthma, RIN, adapter_num, qpcr_conc, bioA_conc, library_tech, num_reads, all.pool.code, orig.pool.code, all.lane.code, all.fc.code)

sum.PC <- prcomp(t(v.gx), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:50] 

##To see if covariates are correlated with a PC (looking at PC1-50)
pval.pca1=matrix(ncol=ncol(gx.covars), nrow=50)
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

setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/CDNull_Processing/")
write.table(pval.pca1, "CDNull_50PCs_vs_covar_TMM_Start_52517.txt", sep="\t", quote=FALSE)

# Modify % variance to show 1 significant digit: 
sumsum.round<-NULL
for(i in 1:length(sumsum$importance[2,])){
  sumsum.round[i]<-round((sumsum$importance[2,i]*100), digits=1)
}

## Plot files: 
# A. CD3 + null, colored by treatment
# B. CD3 + null, colored by genotype: 
pdf("ImHTGX_PCA_CDnullSamp_TMM_RawData_FCreceptor.pdf")
plot(sum.PC$x[,1], sum.PC$x[,2],pch=20, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, main="Colored by FC genotype")   
text(sum.PC$x[,1],sum.PC$x[,2], treatment, cex = .5, pos=3) 
plot(sum.PC$x[,1], sum.PC$x[,2],pch=pheno.file$treat.shape, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 
plot(sum.PC$x[,1], sum.PC$x[,2],pch=pheno.file$treat.shape, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 
text(sum.PC$x[,1],sum.PC$x[,2], indiv, cex = .5, pos=3) 
plot(sum.PC$x[,1], sum.PC$x[,2],pch=20, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, main="Colored by FC genotype")    
text(sum.PC$x[,1],sum.PC$x[,2], fc.geno.letter, cex = .5, pos=3) 
dev.off()


###### 7.  Correct for pool code =========
## Correct for pool effect using removeBatchEffects in limma 
v.nt.nopool<-removeBatchEffect(v.nt, batch=orig.pool.code)

## PCA (no need to change covars)
v.gx<-v.nt.nopool

sum.PC <- prcomp(t(v.gx), scale=FALSE, center=TRUE)
sumsum <- summary(sum.PC)
prop.var<-sumsum$importance
prop.var[1:3,1:50] 

##To see if covariates are correlated with a PC (looking at PC1-50)
pval.pca1=matrix(ncol=ncol(gx.covars), nrow=50)
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

write.table(pval.pca1, "CDNull_50PCs_vs_covar_TMMnoPool.txt", sep="\t", quote=FALSE)
write.table(v.nt.nopool, "VoomTMM_NullCD_Elist_RMPool.txt", sep="\t", quote=FALSE)

# Modify % variance to show 1 significant digit: 
sumsum.round<-NULL
for(i in 1:length(sumsum$importance[2,])){
  sumsum.round[i]<-round((sumsum$importance[2,i]*100), digits=1)
}

## Plot files: 
pdf("ImHTGX_PCA_CDnullSamp_TMM_noPool_FCrecpetor.pdf")
plot(sum.PC$x[,1], sum.PC$x[,2],pch=pheno.file$treat.shape, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 
plot(sum.PC$x[,1], sum.PC$x[,2],pch=pheno.file$treat.shape, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, cex=1.3, main="Colored by FC genotype") 
text(sum.PC$x[,1],sum.PC$x[,2], indiv, cex = .5, pos=3) 
plot(sum.PC$x[,1], sum.PC$x[,2],pch=20, xlab=paste("PC 1 -", (sumsum.round[1]),"% of variance", sep=" "), ylab=paste("PC 2 -",(sumsum.round[2]),"% of variance", sep=" "),col=pheno.file$fc.pca.color, main="Colored by FC genotype")    
text(sum.PC$x[,1],sum.PC$x[,2], fc.geno.letter, cex = .5, pos=3) 
dev.off()
