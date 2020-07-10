#!/bin/bash

set -e 

# Quoted need to be parametrized
# Loop through multiple datasets
# Pass output to rmarkdown and generate pdf
# Document

usage () {
  printf "USAGE:\n"
  printf "`basename $0` <analysis-name> <input-bam-file> <nThreads> <input-hotspots-bed-file (optional)>\n\n"
}

exitWithError () {
  printf "\n***ERROR: $1\n\n"
  usage
  exit
}

if [[ $# == 3 ]]; then
    WORKFLOW="local_workflow.sh"
elif [[ $# == 4 ]]; then
    WORKFLOW="local_workflow_hotspots.sh"
else
    exitWithError "Wrong number of parameters."
fi

export ANALYSIS_NAME="$1"
export INPUT_BAM_name="$2"
export THREADS="$3"
export HOTSPOTS_BED_name="$4"

export MOUNT_DIR="/data"
export LOCAL_DIR="/media/galadriel/fleming/oncopmnet/finalized"
export HGREF="$MOUNT_DIR/ref_genome/hg19/hg19.fasta"
export PARAMS_FILE="$MOUNT_DIR/params/gsfak.json"
export TARGETS_BED="$MOUNT_DIR/bed/solid_custom_merged.bed"
export INPUT_BAM="$MOUNT_DIR/datasets/$ANALYSIS_NAME/$INPUT_BAM_name"
export HOTSPOTS_BED="$MOUNT_DIR/hotspots/$HOTSPOTS_BED_name"

export PREFIX=Sample_$(/media/galadriel/hybridstat/elixir_project/tools/samtools-1.9/samtools view -H "$LOCAL_DIR/datasets/$ANALYSIS_NAME/$INPUT_BAM_name" |fgrep @RG | fgrep SM:| head | tr '\t' '\n' | fgrep SM: | sed -e 's/SM://' -e 's|[ /:$#]|_|g' | uniq)
#export PREFIX=/media/galadriel/hybridstat/elixir_project/tools/samtools-1.9/samtools

# --name Anlysis specific name or ID
# -d detached mode, runs container in the background
# -t allocate a pseudo tty
# -v mounts the volume "finalized" into container's "$MOUNT_DIR" directory
docker run \
--name tmp_"$ANALYSIS_NAME"_"$INPUT_BAM_name" \
-d \
-t \
-v /media/galadriel/fleming/oncopmnet/finalized:"$MOUNT_DIR" sgsfak/tmap-tvc

# -e env variables passed into the container
docker exec \
-e PREFIX \
-e HGREF \
-e PARAMS_FILE \
-e TARGETS_BED \
-e HOTSPOTS_BED \
-e INPUT_BAM_name \
-e ANALYSIS_NAME \
-e INPUT_BAM \
-e THREADS \
-it tmp_"$ANALYSIS_NAME"_"$INPUT_BAM_name" "$MOUNT_DIR"/scripts/"$WORKFLOW"

# Normalize VCF with BCFtools (-norm)
echo Normalizing vcf and breaking multiallelic variants...
"$LOCAL_DIR"/tools/bcftools-1.10.2/bcftools norm \
-f "$LOCAL_DIR"/ref_genome/hg19/hg19.fasta \
-m- "$LOCAL_DIR"/tvc-out/"$ANALYSIS_NAME"/"$PREFIX"/TSVC_variants.vcf \
-o "$LOCAL_DIR"/tvc-out/"$ANALYSIS_NAME"/"$PREFIX"/TSVC_variants_norm.vcf

# Break mutliallelic variants with BCFtools (-norm -m)
# Annotate VCF with Bioconductor's VariantAnnotation
# Add automated hpmnet report generation with Rmarkdown 
# Create relevant repository & documentation

echo "Stopping temporary docker container:"
docker stop tmp_"$ANALYSIS_NAME"_"$INPUT_BAM_name"
echo "Removing temporary docker container:"
docker rm tmp_"$ANALYSIS_NAME"_"$INPUT_BAM_name"
