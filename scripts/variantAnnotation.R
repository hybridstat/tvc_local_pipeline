library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(knitr)
library(stringr)
library(xtable)
library(pander)
library(ggplot2)
library(timeSeries)
library(readr)
library(dplyr)
library(data.table)
library(httr)
library(jsonlite)
library(tidyr)
library(grid)
library(kableExtra)

vcf <- readVcf(vcfDir, "hg19")
### GET REPORT INFO FIELDS ###
sampleId <- "N/A"
analysisId <- "N/A"
sample_name <- sampleName
totalVariants <- length(vcf)

### ANNOTATE VARIANTS ###
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb) <- seqlevels(txdb)[1:24]
annotated <- predictCoding(vcf, txdb, seqSource=Hsapiens)

### CREATE A CUSTOM INPUT (SNV) ###
# Populate SNV input with required values from VCF file
SNV <- data.frame(matrix(ncol = 14, nrow = length(annotated$GENEID)))
colnames(SNV) <- c("Hugo_Symbol", "Variant_Classification", "Protein_Change", "VAR_id", "NCT_IDs", 
"COSMIC ID", "Biomarker", "CDS Mutation", "MAF", "MAC", "Read Depth (DP)", "Genotype", "Allelic Ballance", "Alternate Allele Cov (AO)")
  
aaREF  <- as.vector(annotated$REFAA)
aaVAR  <- as.vector(annotated$VARAA)
aaPOS  <- as.vector(annotated$PROTEINLOC)
cdsREF <- as.character(annotated$REF)
cdsVAR <- as.vector(annotated$varAllele)
cdsPOS <- as.character(annotated$CDSLOC)
  
SNV$Hugo_Symbol <- annotated$GENEID
HGNC_map <- AnnotationDbi::select(org.Hs.eg.db, keys=SNV$Hugo_Symbol, columns="SYMBOL", keytype="ENTREZID")
SNV$Hugo_Symbol <- HGNC_map[match(SNV$Hugo_Symbol,HGNC_map$ENTREZID),2]
  
SNV$Variant_Classification <- as.vector(annotated$CONSEQUENCE)
SNV$Protein_Change <- paste(aaREF ,aaPOS, aaVAR, sep = "")
SNV$VAR_id <- names(annotated)
SNV$Biomarker <- ifelse(SNV$VAR_id == "", 
    paste0(SNV$Hugo_Symbol,":",SNV$Protein_Change," (",SNV$`COSMIC ID`,")",sep = ""), 
    paste0(SNV$Hugo_Symbol,":",SNV$Protein_Change," (",SNV$VAR_id,")",sep = "")
)
SNV$`CDS Mutation` <- paste0("c.", cdsPOS, cdsREF, ">", cdsVAR)
SNV$MAF <- paste(as.vector(info(vcf)[match(SNV$VAR_id,rownames(info(vcf))),13]))
SNV$`Read Depth (DP)` <- (geno(vcf)$DP[match(SNV$VAR_id,rownames(info(vcf)))])
SNV$Genotype <- (geno(vcf)$GT[match(SNV$VAR_id,rownames(info(vcf)))])
SNV$`Alternate Allele Cov (AO)` <- paste(as.vector(info(vcf)[match(SNV$VAR_id,rownames(info(vcf))),5]))
  
SNV <- SNV[!grepl("\\bsynonymous\\b", SNV$Variant_Classification),]
SNV <- unique(SNV)
 
SNV$Variant_Classification <- gsub(pattern = "nonsynonymous", replacement = "missense_variant", SNV$Variant_Classification)
SNV$Variant_Classification <- gsub(pattern = "frameshift", replacement = "frameshift_variant", SNV$Variant_Classification)
SNV$Variant_Classification <- gsub(pattern = "nonsense", replacement = "nonsense_variant", SNV$Variant_Classification)
  
SNV$Protein_Change <- ifelse(test = grepl("frameshift_variant", SNV$Variant_Classification), yes = paste0(SNV$Protein_Change,"FS"), no = SNV$Protein_Change)
