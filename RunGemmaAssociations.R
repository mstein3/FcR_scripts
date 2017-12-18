# This script runs GEMMA in order to identify associations between cytokine levels and rs1801274

library (qvalue)

#reads in arguments from the command line -> Rscript RunGemmaAssociations.R <pheno_list_filename> <pheno column number> <ld_prune_snp_list_file> <plink_file_prefix> <covariate_file> <kinship_file> 
args <- commandArgs(trailingOnly=T)

pheno_list_file=args[1]
#phenotype column number (after the 5th column in fam file -> i=1 is column 6)
i=as.numeric(args[2])
ld_prune_snp_list_file=args[3]
plink_file_prefix=args[4]
covariate_file=args[5]
kinship_file=args[6]


#read in pheno list to label outputfile
pheno_list=read.table(pheno_list_file,sep="\t",header=F)
phenos=pheno_list[,1]

#set output directory 
setwd("/path/to/output/dir")

#run_gemma from R
system(paste("gemma -bfile ", plink_file_prefix, " -k ",kinship_file ," -c ", covariate_file, " -n ", i+1, " -km 2 -lmm 1 -notsnp -miss 1 -o ", phenos[i], sep=""))





