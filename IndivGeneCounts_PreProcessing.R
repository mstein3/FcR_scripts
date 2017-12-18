#### Processing RNAseq gene counts after mapping with STAR

####### 1. Count the number of mapped reads =============
setwd("~/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/mapped_reads_outputs_hg19/")
file_list <- list.files()
map_read_counts<-data.frame(matrix(NA, nrow = length(file_list), ncol = 5))

for(i in 1:length(file_list)) {
    dataset <- read.table(file_list[i], sep="\t")
    total.reads<-sum(dataset[,2])
    unmap.reads<-sum(dataset[1:4,2])
    map.reads<- (total.reads-unmap.reads)
    percent.map<- (map.reads/total.reads)
    unmap.reads.strict<-dataset[1,2]
    map.reads.strict<- (total.reads-unmap.reads.strict)
    percent.map.strict<-(map.reads.strict/total.reads)
    map_read_counts[i,1]<-paste(file_list[i])
    map_read_counts[i,2]<-map.reads
    map_read_counts[i,3]<-percent.map
    map_read_counts[i,4]<-map.reads.strict
    map_read_counts[i,5]<-percent.map.strict
  
}

map.table.names<-c("name", "NumMapReads_allAmbig", "PercentMapReads_allAmbig", "NumMapReads_UnmapOnly", "PercentMapReads_UnmapOnly")
colnames(map_read_counts)<-map.table.names

setwd("~/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/")
write.table(map_read_counts, "Mapped_read_counts.txt", sep="\t", quote=FALSE, row.names=FALSE)

####### 2. Read in the gene count files for processing. ========
setwd("~/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/mapped_reads_outputs_hg19/")

file_list_test<-file_list[1:5]

for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, sep="\t")
  }
  
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-read.table(file, sep="\t")
    colnames(temp_dataset)<-c("V1", paste(file), "V3", "V4")
    dataset<-cbind.data.frame(dataset, temp_dataset[,2])
    rm(temp_dataset)
  }
  
}

dim(dataset)
# [1] 57824  1372

### Remove columns 2, 3 and 4 (from the first file read in)
dataset2<-dataset[,-2:-4]

# Read in file list with unique extract ord names, all with 4 digits. 
# Name columns by cleaning up some of the file list: 
file_list4<-read.table("file_list4.txt", sep="\t", header=FALSE)
colnames(dataset2)<-file_list4$V1

# Write this file 
setwd("~/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R/")
write.table(dataset2, "GeneCounts_SamplesNotMerged.txt", sep="\t", row.names=F)

# Remove the first 4 rows from the dataset, which do not correspond to genes:
dataset4<-dataset2[-1:-4,]

# Total up number of reads for each seq run: 
each.run.map<-colSums(dataset4)

###### 3.  Add together the technical replicates =============
colnames_for_org <- unique(sapply(strsplit(names(dataset4), "_", fixed=TRUE), "[", 1) )
# Add together columns with matching extract ords, using the uniform 4-digit codes (ex: 0007, 0087, etc.)
res<-sapply(colnames_for_org, function(x) rowSums(dataset4 [, grep(paste(x,"_",sep=""), names(dataset4)), drop=FALSE] )  )

# Write out combined gene count table
setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R")
write.table(res, "RawGeneCounts_RunsCombined_ENSG_names_noQC.txt", sep="\t", quote=FALSE)

####### 4. Remove X, Y, M genes ==================
setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing/Processing_in_R")
ensg.all<-read.table("ENSG_GeneNames_Chr_All.txt", sep=" ", header=FALSE)

# Subset out protein coding genes:
protcode<-ensg.all$V5=="protein_coding"
ensg.pc<-ensg.all[protcode,]
dim(ensg.pc)
#[1] 20345     6

### 4.1 Remove X, Y, and M chr: 
xym.chr<-ensg.pc$V1=="chrX" | ensg.pc$V1=="chrY" | ensg.pc$V1=="chrM"
ensg.auto<-ensg.pc[!xym.chr,]

# Subset: Subset out protein coding + autosomal genes + genes from res from ensg file: 
ensg.intersect<-intersect(rownames(res), ensg.auto[,4])
ensg.match<-match(ensg.intersect, ensg.auto[,4])
ensg.auto2<-ensg.auto[ensg.match,]

# Subset: Subset out protein coding + autosomal genes + genes from res from res file: 
res.match<-match(ensg.intersect, rownames(res))
res2<-res[res.match,]

######### 5. Remove samples with <7M mapped reads: =========
setwd("/Users/michellestein/Dropbox/Winter_2012_studies/Full_Hutterite/ImHT_GX/2017_RNAseq_processing")
pheno.start<-read.table("Pheno_Rready_403comb.txt", header=TRUE, sep="\t")
enough.reads<-pheno.start$Total_Map>7000000
res3<-res2[,enough.reads]
pheno_7M<-pheno.start[enough.reads,]

######## 6. Split up files by treatment for downstream processing: ===============
null1<-pheno_7M$treat=="null"
cd1<-pheno_7M$treat=="CD"

null.pheno<-pheno_7M[null1,]
cd.pheno<-pheno_7M[cd1,]

null.gc<-res3[,null1]
cd.gc<-res3[,cd1]

# Combine into two sets of treatments: 
nt.pheno<-rbind.data.frame(null.pheno, cd.pheno)
nt.gc<-cbind.data.frame(null.gc, cd.gc)

# Write out all files 