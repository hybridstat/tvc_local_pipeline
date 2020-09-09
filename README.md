# TVC Local Pipeline

NOT TO BE MERGED WITH MASTER.
This is a branch that offers the option to generate a report given the three oncoPMNET report input CSVs

## Usage

```sh
$ ./oncompnet.sh <analysis-name> <input-bam-file> <nThreads> <input-hotspots-bed-file (optional)>
```

## Required Tools

The following tools are required

* BCFTools - Normalize vcf and breake multiallelic variants
* Docker - Container which runs the main part of the analysis
* R version 3.6.2 - Variant annotation and report generation
    * VariantAnnotation
    * BSgenome.Hsapiens.UCSC.hg19
    * TxDb.Hsapiens.UCSC.hg19.knownGene
    * org.Hs.eg.db
    * knitr
    * stringr
    * xtable
    * pander
    * ggplot2
    * timeSeries
    * readr
    * dplyr
    * data.table
    * httr
    * jsonlite
    * tidyr
    * grid
    * kableExtra
* WeasyPrint - Convert HTML rmarkdown output to pdf report

## Directory Format

| Directory | Content |
| ------ | ------ |
| **datasets** | Contains the analyses subdirs. Each analysis dir can contain multiple BAM files belonging to each analysis|
| **params** |  TVC parameters json file presets. Currently hard-coded in oncopmnet shell script |
| **bed** |  PMNet custom BED file |
| **hotspots** | Hotspot BED files to be used with the hotspot option |
| **ref_genome** | Hg19 reference genome files |
| **tvc-out** | Contains the output files organized in the same way as in datasets dir |
| **scripts** | Scripts used in the analysis |

## Running an analysis
#### Dataset creation
Analysis usually begins from a BAM file which is downloaded from local ION torrent server. The BAM file needs to be copied inside the analysis sub-directory under the datasets dir. Any required analysis sub-dirs will need to be manually created before the analysis. If the BAM file belongs to an already existing analysis sub-dir then it just needs to be copied to the right place.

*WARNING: If the BAM file input in an analysis sub-dir already exists the analysis will exit with error.*
Generally, there can be no two identically named BAM files under the same analysis. If this is the case: 
* either a different name must be given to the BAM file, or
* the BAM file should be assigned to a different analysis

#### Hotspot parameter
The option to provide a hotspots BED file is provided. In this case the name of the hotspot file should be provided as the 5th positional arguement, which is optional. 

*WARNING: The hotspots file must be inside the hotspots directory.*

## Analysis steps
The pipeline follows the next steps to produce the final report:
1. Creation of a temporary container based on this image: https://hub.docker.com/r/sgsfak/tmap-tvc
2. Inside the container
    - Preparation of Hotspots VCF with tvcutils (if needed)
    - Mapping with TMAP
    - Aligned BAM sorting with samtools
    - Sorted BAM indexing with samtools
    - Alignment stats generation with samtools
    - Variant calling with TVC
3. Normalization of vcf and splitting of multi-allelic variants with BCFtools
4. Variant annotation with Bioconductor VariantAnnotation
5. Creation of HTML report with Rmarkdown
6. Convertion of HTML report to PDF format with Weasyprint

## Example
Running an analysis will require the following detailed commands after having followed instructions described in Dataset Creation  chapter:

```sh
$ cd /media/galadriel/fleming/oncopmnet/oncopmnet_pipeline
$ ./oncompnet.sh test_analysis PMN_S1-Q1DP 12 ocp_hotspot_blist_generated.bed>
```
