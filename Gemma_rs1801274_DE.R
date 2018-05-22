#!/usr/bin/env Rscript

# Run GEMMA for differential expression by rs1801274 genotype, using PLINK file format. 

# Usage: Rscript Gemma_rs1801274_DE.R <counts> <covar_file> <relat> <prefix_plink> <gene_name> <dir_pheno> <dir_plink> <dir_gemma> 

# counts: gene-by-sample matrix of normalzied gene expression values
#
# covar: tab-separated file with  covariates.
# First column is 1 for intercept. No column names.
#
# relat: pairwise relatedness values in GEMMA format for the Hutterites
#
# prefix_plink: The prefix of the PLINK files (.bim, .bed, .fam)
#
# gene_name: Name of the gene 
#
# dir_pheno: to save subsetted counts file
#
# dir_plink: to save plink files for each gene
#
# dir_gemma: to save gemma outputs.


# ----------------Read in inputs-----------------------------#
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 8)

# counts
f_exp <- args[1]
# covar
f_covar <- args[2]
# relat
f_relat <- args[3]
# prefix_plink
plink_prefix_all <- args[4]
# gene_name
gene_name<-args[5]
# dir_pheno
dir_pheno <- args[6]
# dir_plink
dir_plink <- args[7]
# dir_gemma
dir_gemma <- args[8]

# Add file types to plink prefix:
f_plink_bed <- paste0(plink_prefix_all, ".bed")
f_plink_bim <- paste0(plink_prefix_all, ".bim")
f_plink_fam <- paste0(plink_prefix_all, ".fam")

stopifnot(file.exists(f_exp, f_covar, f_relat,
                      f_plink_bed, f_plink_bim, f_plink_fam))
                      
# Import data files
exp <- read.delim(f_exp, check.names = FALSE)
covar <- read.table(f_covar)
plink_bim <- read.table(f_plink_bim, stringsAsFactors = FALSE)
plink_fam <- read.table(f_plink_fam, stringsAsFactors = FALSE)

stopifnot(ncol(exp) == nrow(covar))

colnames(covar)<-c("intercept", "age", "sex")

# Make output subdirectories in a separate file.
fc_3geno_subdir <- paste0("fc_3geno_subdir", "/")
dir_pheno_covar <- file.path(dir_pheno, fc_3geno_subdir)
dir_plink_covar <- file.path(dir_plink, fc_3geno_subdir)
dir_gemma_covar <- file.path(dir_gemma, fc_3geno_subdir)

# Create global output file
f_top <- file.path(dir_gemma, paste0("top-fc-geno.txt"))
f_top_colnames <- c("gene", "chr", "rs", "ps", "n_miss", "allele1", "allele0",
                      "af", "beta", "se", "log1_H1", "l_remle", "l_mle", "p_wald", "p_lrt", "p_score", "n_snps")
cat(f_top_colnames, file = f_top, sep = "\t")
cat("\n", file = f_top, sep = "", append = TRUE)  


#----------------- Functions for differential expression with a genotype analysis ---------------------------#

# Phenotype file: col1=FID, col2=IID, col3=gene expression/phenotype
create_pheno_file <- function(fid, iid, e, f = "") {
  stopifnot(length(fid) == 1 | length(fid) == length(e),
            length(iid) == 1 | length(iid) == length(e),
            is.numeric(e))
  d <- data.frame(fid, iid, e)
  write.table(d, file = f, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)
  return(invisible(f))
}

create_plink_g <- function(prefix_in, f_pheno, prefix_out) {
  cmd <- sprintf("plink --bfile %s --make-bed --pheno %s --out %s",
                 prefix_in, f_pheno, prefix_out)
  system(cmd)
  return(invisible(cmd))
}


run_gemma <- function(plink, relatedness, covar, out_prefix, outdir) {
  cmd <- sprintf("gemma -bfile %s -k %s -km 2 -maf 0 -miss 0.2 -c %s -lmm 4 -o %s",
                 plink, relatedness, covar, out_prefix)
  system(cmd)
  # Move to output directory
  cmd2 <- sprintf("mv output/%s* %s", out_prefix, outdir)
  system(cmd2)
  outfile <- paste0(outdir, out_prefix, ".assoc.txt")
  return(invisible(outfile))
}

parse_gemma <- function(f, gene, outfile = "") {
  gemma <- read.delim(f, stringsAsFactors = FALSE)
  n_snps <- nrow(gemma)
  top <- data.frame(gene, gemma[which.min(gemma$p_lrt), ], n_snps,
                    stringsAsFactors = FALSE)
  write.table(top, file = outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  return(invisible(outfile))
}

# Adapted from https://github.com/jdblischak/cardioqtl/blob/master/scripts/run-gemma.R
# ------------------ DE analysis ---------------------------------------#

g <- gene_name
gene_row_num<-which(rownames(exp) == g)
iid <- colnames(exp)
fid <- "HUTTERITES"

message("\n\n####\ngene: ", g, "\n####\n\n")
f_pheno_g <- paste0(dir_pheno_covar, g, ".pheno")

create_pheno_file(fid = fid, iid = iid, e = as.numeric(exp[gene_row_num, ]),
                    f = f_pheno_g)
                    
plink_prefix_g <- paste0(dir_plink_covar, g)

message("Running PLINK")
create_plink_g(prefix_in = plink_prefix_all, f_pheno = f_pheno_g, prefix_out = plink_prefix_g)  

if (!file.exists(paste0(plink_prefix_g, ".bed"))) 
    message("\nskipped:\tNo PLINK output\n") 
                 
message("Running GEMMA")
  f_gemma <- run_gemma(plink = plink_prefix_g, relatedness = f_relat, covar = f_covar,
                       out_prefix = paste0("fc_cd_diff-", g),
                       outdir = dir_gemma_covar)
                       
message("Parse GEMMA")
f_gemma_top <- sub(".assoc.txt", ".top.txt", f_gemma)
parse_gemma(f = f_gemma, gene = g, outfile = f_gemma_top)
# Write to global results file
system(sprintf("cat %s | sed -e '1d' >> %s", f_gemma_top, f_top))
          




