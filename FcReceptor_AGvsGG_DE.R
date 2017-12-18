### Run limma on AG vs GG in CD3-null response: 

### Read in raw gene expression data: 
setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/PreQC_GeneCount_pheno_files/")
nt.pheno<-read.delim("Pheno_Rready_nullCD_7M.txt", sep="\t")
setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/CDNull_Processing/")
nt.gc3<-read.table("RawGC_DetectGenes_NullCD_GeneNames.txt", sep="\t")


##### 1. Remove individuals with AA genotypes =============
aa.genos<-nt.pheno3$FcConsensus_letter=="AA"
nt.pheno4<-nt.pheno3[!aa.genos, ]
nt.gc4<-nt.gc3[,!aa.genos]
dim(nt.pheno4)
#[1] 205  87

####### 2. Find matches the the null and CD3 samples ==============
n1<-nt.pheno4$treat=="null"
pheno.null<-nt.pheno4[n1,]
null.gc<-nt.gc4[,n1]

t1<-nt.pheno4$treat=="CD"
pheno.cd<-nt.pheno4[t1,]
cd.gc<-nt.gc4[,t1]

# Find overlap of findivs between null.gc and cd,gc: 
overlap.5<-intersect(pheno.null$FINDIV, pheno.cd$FINDIV)
overlap.7n<-match(overlap.5, pheno.null$FINDIV)
overlap.7t<-match(overlap.5, pheno.cd$FINDIV)

pheno.null.match<-pheno.null[overlap.7n,]
null.gc.match<-null.gc[,overlap.7n]
pheno.cd.match<-pheno.cd[overlap.7t,]
cd.gc.match<-cd.gc[,overlap.7t]

# Merge these files: 
pheno.ntm<-rbind.data.frame(pheno.null.match, pheno.cd.match)
ntm.gc<-cbind.data.frame(null.gc.match, cd.gc.match)
# Drop subsetted factors: 
pheno.ntm$treat<-factor(pheno.ntm$treat)
pheno.ntm$FcConsensus_let<-factor(pheno.ntm$FcConsensus_let)

####### 3. Run Limma ====================
# Model data using the group means parameterization. Each beta represents the mean expresion level for samples belonging to that group. 
library(limma)
table(pheno.ntm$FcConsensus_let, pheno.ntm$treat)
#      CD null
#  AG 42   42
#  GG 47   47

test_factor<-paste(pheno.ntm$FcConsensus_let, pheno.ntm$treat, sep=".")
test_factor<-factor(test_factor, levels=c("AG.null", "AG.CD", "GG.null", "GG.CD"))

design<-model.matrix(~0 + test_factor + as.factor(pheno.ntm$sex) + as.numeric(pheno.ntm$Age) + as.factor(pheno.ntm$Pool_1.0), groups=pheno.ntm)
colnames(design)<-c("AG.null","AG.CD","GG.null","GG.CD", "sex", "age", "pool1", "pool2", "pool3", "pool4", "pool5", "pool6", "pool7", "pool8", "pool9", "pool10", "pool11", "pool12", "pool13", "pool14", "pool15", "pool16", "pool17", "pool18", "pool19", "pool20", "pool21", "pool22", "pool24")

library(limma)
library(edgeR)

#### Normalize the data 
dge.ntm <- DGEList(counts=ntm.gc)
dge.ntm <- calcNormFactors(dge.ntm)

# Apply the voom transformation 
v.ntm <- voom(dge.ntm, design=design,plot=TRUE)

# Create the individual blocks
indiv_effect <- duplicateCorrelation(v.ntm, design, block = pheno.ntm$FINDIV)

fit <- lmFit(v.ntm, design, block = pheno.ntm$FINDIV,
             correlation = indiv_effect$consensus.correlation)
cont_mat <- makeContrasts(AG.resp = AG.CD - AG.null,
                          GG.resp = GG.CD - GG.null,
                          diff = (GG.CD - GG.null) -
                            (AG.CD - AG.null),
                          levels = design)
fit2 <- contrasts.fit(fit, cont_mat)
fit2 <- eBayes(fit2)
bh_fdr_cutoff <- 0.05
results <- decideTests(fit2, p.value = bh_fdr_cutoff)
pdf("AG_GG_CD3resp_limma.pdf")
vennDiagram(results, main="FDR 5%")
dev.off()

# Save this file 
saveRDS(fit2, "Limma_array_CD3resp_Fc.rds")
fit.toptable<-topTable(fit2,sort="none",n=Inf)
write.table(fit.toptable, "TopTable_allgenes_CD3resp_Fc.txt", sep="\t", quote=FALSE)

# Pull out P-values and coefficients, and contrast matrix:   
head(fit2$p.value)

all.pvals<-fit2$p.value
all.coef<-fit2$coefficients

###### 4. Separate out DE genes =======================
## Separate out different groups of genes

# A. Genes that respond in AG, GG, and respond differently
diff.resp<-which(results[,1]!=0 & results[,2]!=0 & results[,3]!=0)
length(diff.resp)
# [1] 2197

# B. Genes that respond similarly in both AG and GG 
same.resp<-which(results[,1]!=0 & results[,2]!=0 & results[,3]==0)
length(same.resp)
# [1] 6589

# C. Genes that respond in AG, but not GG, not differently. 
AG.resp<-which(results[,1]!=0 & results[,2]==0 & results[,3]==0)
length(AG.resp)
# [1] 265

# D. Genes that respond in GG, but not AG, not differently. 
GG.resp<-which(results[,1]==0 & results[,2]!=0 & results[,3]==0)
length(GG.resp)
#[1] 411
 
# Genes that respond differently, and not to AG or GG: 
one.diff<-which(results[,1]==0 & results[,2]==0 & results[,3]!=0)
# AIG1 

# Genes that respond to AG, and differently: 
AG.diff<-which(results[,1]!=0 & results[,2]==0 & results[,3]!=0)
length(AG.diff)
# [1] 41

# Genes that respond to GG, and differently: 
GG.diff<-which(results[,1]==0 & results[,2]!=0 & results[,3]!=0)
length(GG.diff)
# [1] 91

# For the groups with small numbers of genes, subset from fit.toptable, and write to file: 
# Genes that respond differently, and not to AG or GG: 
one.diff.tab<-fit.toptable[one.diff,]
gg.diff.tab<-fit.toptable[GG.diff,]
ag.diff.tab<-fit.toptable[AG.diff,]

# For the groups with more genes, write out to file. 
gg.resp.tab<-fit.toptable[GG.resp,]
ag.resp.tab<-fit.toptable[AG.resp,]
same.resp.tab<-fit.toptable[same.resp,]
diff.resp.tab<-fit.toptable[diff.resp,]
# Write out these files 

# Break down the diff response genes: 
# Where AG is 1, and GG is 1
pos.resp<-which(results[,1]==1 & results[,2]==1 & results[,3]!=0)
length(pos.resp)
# [1] 1393

# Where AG is -1 and GG is 1 
neg.resp<-which(results[,1]==-1 & results[,2]==-1 & results[,3]!=0)
length(neg.resp)
#[1] 801

rev.resp<-which(results[,1]==-1 & results[,2]==1 & results[,3]!=0)
#PAQR7 MIER2 

# Where is the magnitude of the response greater in GG? 
gg.mag<-which(abs(diff.resp.tab$GG.resp)>=abs(diff.resp.tab$AG.resp))
length(gg.mag)
# [1] 2150

## Pull out the 2197 genes that respond differently in AG and GG: 
diff.resp.v<-v.ntm[which(results[,1]!=0 & results[,2]!=0 & results[,3]!=0),]

####### 5. Volcano plots for CD3 response in GG. ================
# Read in colors: 
dgreen<-c("#1B9E77") # GG
#rgb: 27, 158, 119
dorange<-c("#D95F02") # AA
# rgb: 271, 95, 2
dpurple<-c("#7570B3") #AG
# rgb: 117, 112, 179 

#dgreengray<-c("#4E6A1") ## GG, up in Null
#dgreengray<-rgb(92, 92, 92, maxColorValue = 255)
# rgb: 78, 106, 97

#dpurplegray<-c("#88869c") # AG, up in Null
# rgb: 136, 134, 156
#dvolcpurple<-c("#5d55ce") # AG, up in CD. Higher saturation than std color. 
# rgb: 93, 85, 206 

### Set up elements needed for volcano plots: 
library(ggplot2)

## Split up samples to get means 
gx.ntm<-v.ntm$E

# Separate out by genotype: 
gg1<-pheno.ntm$FcConsensus_letter=="GG"
pheno.mg<-pheno.ntm[gg1,]
gx.mg<-gx.ntm[,gg1]

# Separate out by treatment: 
null1<-pheno.mg$treat=="null"
pheno.mgn<-pheno.mg[null1,]
pheno.mgt<-pheno.mg[!null1,]
gx.mgn<-gx.mg[,null1]
gx.mgt<-gx.mg[,!null1]

# Make x-axis values: 
mean_mgn<-apply(gx.mgn, 1, mean)
mean_mgt<-apply(gx.mgt, 1, mean)
x.axis.gg<- mean_mgt- mean_mgn
fc.gg<-2^x.axis.gg

y.axis.gg<- (-log10(all.pvals[,2]))
xy.axis.gg<-cbind.data.frame(x.axis.gg, y.axis.gg, names(y.axis.gg))
colnames(xy.axis.gg)<-c("x.axis.gg", "y.axis.gg", "Gene")

# Use results parameter to find genes that are DE. 
xy.axis.gg$color_flag<-ifelse(results[,2]=="0", "black", ifelse(results[,2]=="1", dgreen, "gray50"))

gg<-ggplot(data=xy.axis.gg, aes(x=x.axis.gg, y=y.axis.gg, color=color_flag)) +
  geom_point(shape=16, colour=xy.axis.gg$color_flag, size=1.8)+
  ylab("-log10(P value)")+
  xlab("log2(fold change)")+
  geom_vline(xintercept=0, color="gray15")+
  ylim(0,128) +
  xlim(-8, 8) +
  ggtitle("GG")+
  #theme_bw()+
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), axis.title.x=element_text(size=20), axis.title.y=element_text(size=18))
gg

# Install ggrepel
library(ggrepel)
gg + geom_text_repel(data=subset(xy.axis.gg, x.axis.gg > 5 | x.axis.gg < -5), aes(label=Gene), colour="black", size=4)

pdf("CD3_resp_GG_VolcPlot.pdf")
gg + geom_text_repel(data=subset(xy.axis.gg, x.axis.gg > 5 | x.axis.gg < -5), aes(label=Gene), colour="black", size=4)
dev.off()

########## 6 Make Volcano plot for AG response to CD3 treatent. ========
# Separate out by genotype: 
ag1<-pheno.ntm$FcConsensus_letter=="AG"
pheno.ma<-pheno.ntm[ag1,]
gx.ma<-gx.ntm[,ag1]

# Separate out by treatment: 
null2<-pheno.ma$treat=="null"
pheno.man<-pheno.ma[null2,]
pheno.mat<-pheno.ma[!null2,]
gx.man<-gx.ma[,null2]
gx.mat<-gx.ma[,!null2]

# Make x-axis values: 
mean_man<-apply(gx.man, 1, mean)
mean_mat<-apply(gx.mat, 1, mean)
x.axis.ag<- mean_mat- mean_man
fc.ag<-2^x.axis.ag

y.axis.ag<- (-log10(all.pvals[,1]))
xy.axis.ag<-cbind.data.frame(x.axis.ag, y.axis.ag, names(y.axis.ag))
colnames(xy.axis.ag)<-c("x.axis.ag", "y.axis.ag", "Gene")

xy.axis.ag$color_flag<-ifelse(results[,1]=="0", "black", ifelse(results[,1]=="1", dpurple, "gray50"))

ag<-ggplot(data=xy.axis.ag, aes(x=x.axis.ag, y=y.axis.ag, color=color_flag)) +
  geom_point(shape=16, colour=xy.axis.ag$color_flag, size=1.8)+
  ylab("-log10(P value)")+
  xlab("log2(fold change)")+
  geom_vline(xintercept=0, color="gray15")+
  ylim(0,128) +
  xlim(-8, 8) +
  ggtitle("AG")+
  #theme_bw()+
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), axis.title.x=element_text(size=20), axis.title.y=element_text(size=18))
ag
ag + geom_text_repel(data=subset(xy.axis.ag, x.axis.ag > 5 | x.axis.ag < -5), aes(label=Gene), colour="black", size=4)

pdf("CD3_resp_AG_VolcPlot.pdf")
ag + geom_text_repel(data=subset(xy.axis.ag, x.axis.ag > 5 | x.axis.ag < -5), aes(label=Gene), colour="black", size=4)
dev.off()

####### 7. Make a Beta value difference plot ===============
diff.genes.all<-c(rownames(diff.resp.tab), rownames(gg.diff.tab), rownames(ag.diff.tab))
diff.genes.all<-c(rownames(diff.resp.tab))
length(diff.genes.all)

# Combine x.axis.gg and x.axis.ag
fc.diff.all<-cbind.data.frame(x.axis.gg, x.axis.ag)
match.diff<-match(diff.genes.all, rownames(fc.diff.all))
fc.diff<-fc.diff.all[match.diff,]
colnames(fc.diff)<-c("log2.fc.gg", "log2.fc.ag")

## To plot avg gene expression, use the log2 fold change against each other in just the genes of interest: 
plot(fc.diff[,1], fc.diff[,2], pch=19, xlab="gg", ylab="ag", col="blue")
abline(h=0, col="black")
abline(v=0, col="black")
abline(a=0, b=1, col="gray60")

fc.diff$color_flag<-ifelse(abs(fc.diff$log2.fc.gg)>abs(fc.diff$log2.fc.ag), dgreen, dpurple)

fc<-ggplot(data=fc.diff, aes(x=log2.fc.gg, y=log2.fc.ag, color=color_flag)) +
  geom_point(shape=16, colour=fc.diff$color_flag, size=1.8)+
  ylab("AG \n log2(fold change)")+
  xlab("GG \n log2(fold change)")+
  geom_vline(xintercept=0, color="gray10")+
  geom_abline(slope=1, intercept=0, colour="gray50")+
  geom_hline(yintercept=0, color="gray10")+
  ylim(-8, 8) +
  xlim(-8, 8) +
  ggtitle("Gene expression fold change responses")+
  #theme_bw()+
  theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16), axis.title.x=element_text(size=18), axis.title.y=element_text(size=16))
fc

pdf("CD3_resp_CompareGXResp2197.pdf")
fc
dev.off()

## Count the number of genes in each quadrant: 
# Positive and increased in GG: 
pos.gg<-fc.diff$color_flag=="#1B9E77" & fc.diff$log2.fc.gg > 0
# 1397 genes
# with just 2197, 1357 genes
# Negative and larger magnitude in GG: 
neg.gg<-fc.diff$color_flag=="#1B9E77" & fc.diff$log2.fc.gg < 0
# 844 genes
# with just 2197, 793

# Positive and increased in AG: 
pos.ag<-fc.diff$color_flag!="#1B9E77" & fc.diff$log2.fc.ag > 0
# 60 genes
# with just 2197, 39 genes
neg.ag<-fc.diff$color_flag!="#1B9E77" & fc.diff$log2.fc.ag < 0
# 28 genes 
# with just 2197, 8 