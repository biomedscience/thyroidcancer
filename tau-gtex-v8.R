# BiocManager::install("biomaRt")
BiocManager::install("EnsDb.Hsapiens.v86")

library("openxlsx")
library('biomaRt')
library("AnnotationHub")
library("EnsDb.Hsapiens.v86")
library("tidyverse")

data<-read.csv("./extdata/GSI/Tau_gene_V8.csv")
head(data)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- data$gene_id
G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=genes,mart= mart)

ensembl.genes <- data$gene_id
geneID <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))

input<-left_join(geneID,data,by=c("GENEID"="gene_id"))
#input %>% dplyr::filter(SYMBOL=="AIM2")

metainput<-read.table(file=".//extdata//GSE55457_GSE55584_GSE55235_RA_GPL960_meta.txt")
metainput <-data.table(metainput,keep.rownames = TRUE)
metainput<-metainput %>% dplyr::filter(pval<0.0001)

output<-left_join(metainput,input,by=c("rn"="SYMBOL"))
head(output)
output %>% dplyr::filter(rn=="AIM2")
output<-output %>% dplyr::filter(tau>0.7,pval<0.000001)
write.csv(output,file="meta.tau.csv")
