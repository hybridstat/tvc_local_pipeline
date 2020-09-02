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

### MATCH CANCER TYPES ###
# Match cancer type between input and available synonyms
synonyms = read.csv("/media/galadriel/fleming/oncopmnet/oncopmnet_pipeline/report_extra/cancer_types.csv", header = TRUE,sep=",")

if (!exists("cancer") || !any(cancer == synonyms$tcga_cancer)) {
  warning(paste0('Cancer acronym not found. Using unspecified cancer type as input'))
  cancer <- "N/A"
}
cancer_report    = unique(as.character(synonyms[grep(paste0("\\b", cancer, "\\b"),synonyms$tcga_cancer,ignore.case = T),"Report"]))


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
SNV <- data.frame(matrix(ncol = 7, nrow = length(annotated$GENEID)))
colnames(SNV) <- c("Γονίδιο", "Γεννετική Παραλλαγή", "Κλινική σημασία", "Συχνότητα παραλλαγμένου αλληλόμορφου (VAF)", "Συνολική κάλυψη θέσης (coverage)", "Variant_Classification", "Protein_Change")
  
aaREF  <- as.vector(annotated$REFAA)
aaVAR  <- as.vector(annotated$VARAA)
aaPOS  <- as.vector(annotated$PROTEINLOC)
cdsREF <- as.character(annotated$REF)
cdsVAR <- as.vector(annotated$varAllele)
cdsPOS <- as.character(annotated$CDSLOC)
  
SNV[,1] <- annotated$GENEID
HGNC_map <- AnnotationDbi::select(org.Hs.eg.db, keys=SNV[,1], columns="SYMBOL", keytype="ENTREZID")
SNV[,1] <- HGNC_map[match(SNV[,1],HGNC_map$ENTREZID),2]
  
SNV[,6] <- as.vector(annotated$CONSEQUENCE)
SNV[,7] <- paste(aaREF ,aaPOS, aaVAR, sep = "")
SNV$VAR_id <- names(annotated)
SNV[,2] <- paste0(SNV[,1],":",SNV$Protein_Change,sep = "")

SNV[,3] <- "Pathogenic"

SNV[,4] <- paste(as.vector(info(vcf)[match(SNV$VAR_id,rownames(info(vcf))),13]))
SNV[,5] <- paste(as.vector(info(vcf)[match(SNV$VAR_id,rownames(info(vcf))),5]))
  
SNV <- SNV[!grepl("\\bsynonymous\\b", SNV$Variant_Classification),]
SNV$Variant_Classification <- gsub(pattern = "nonsynonymous", replacement = "missense_variant", SNV$Variant_Classification)
SNV$Variant_Classification <- gsub(pattern = "frameshift", replacement = "frameshift_variant", SNV$Variant_Classification)
SNV$Variant_Classification <- gsub(pattern = "nonsense", replacement = "nonsense_variant", SNV$Variant_Classification)
SNV$Protein_Change <- ifelse(test = grepl("frameshift_variant", SNV$Variant_Classification), yes = paste0(SNV$Protein_Change,"FS"), no = SNV$Protein_Change)


annotable <- data.frame(matrix(ncol = 2, nrow = nrow(SNV)))
colnames(annotable) <- c("Γενωμική παραλλαγή", "Ταξινόμηση Παραλλαγής")

annotable[,1] <- ifelse(SNV$VAR_id == "", 
    paste0(SNV[,1],":",SNV$Protein_Change," (",SNV$`COSMIC ID`,")",sep = ""), 
    paste0(SNV[,1],":",SNV$Protein_Change," (",SNV$VAR_id,")",sep = "")
)
annotable[,2] <- SNV[,6]


metrics_file <- list.files(pattern='alignment_stats.txt', path=output, full.names=T)
metrics <- read.delim(file=metrics_file, sep="\t", header=F)
colnames(metrics) <- c("Συνολικά Διαβάσματα", "Στοιχισμένα Διαβάσματα", "Στοιχισμένα Διαβάσματα (>Q20)", "Μέσο μήκος διαβάσματος")

results <- unique(SNV[,1:5])
rownames(results) <- NULL

annotations <- unique(annotable)
rownames(annotations) <- NULL

signatures <- data.frame(matrix(ncol = 5, nrow = 1))
colnames(signatures) <- NULL
rownames(signatures) <- NULL
signatures[,] <- ""
signatures[1,1] <- "Υπεύθυνος Ανάλυσης"
signatures[1,5] <- "Επιστημονικός Υπεύθυνος"
