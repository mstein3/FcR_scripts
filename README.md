# FcR_scripts
Scripts used to investigate the effect of rs1801274 on anti-CD3 + anti-CD28 treated whole blood

## 1. Cytokine analysis: 
RunAssociationScript.sh: runs the .pbs script for all cytokines in the study, simultaneously. 

StartGemmaRun.pbs: Submits the .R script with the files required (phenotype file, covariate file, kinship file, etc). 

RunGemmaAssociations.R: Runs GEMMA (http://stephenslab.uchicago.edu/software.html) to calculate associations between rs1801274 and cytokine levels. GEMMA acocunts for relatedness with a genotype-based kinship matrix. 

## 2. Gene Expression, data processing: 
fastqc_run.pbs: Runs fastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to QC RNAseq data.

multiQC_run.pbs: Runs multiQC (http://multiqc.info/), a tool to aggreate and visualize fastQC outputs. 

verifyBAM.pbs: Runs verifyBamID (https://genome.sph.umich.edu/wiki/VerifyBamID) to identify possible sample swaps and/or contamination when external genotype data is available. 

star_run.pbs: Runs STAR(https://github.com/alexdobin/STAR), to align reads to hg19 genome, and counts genes for each sample. 

IndivGeneCounts_PreProcessing.R: Organize gene counts, add together sequencing runs from the same sample, remove genes that map to chr X, Y, and mitochondrial genes, and remove non-protein-coding genes. 

GeneCounts_Processing.R: Further processing to remove genes with low count per million, principal components analysis (PCA) to identify confounding technical covariates.

## 3. Gene Expression, data analysis: 
FcReceptor_GX_PCA.R: Plot gene expression data using PCA.

FcReceptor_AGvsGG_DE.R: Calculate differential expression of genes between individuals with genotype AG and individuals with genotype GG for rs1801274.

