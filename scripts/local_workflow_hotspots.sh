#!/bin/bash

printToLog() {
	date +"%F__%T ---> $1"
}

usage () {
  printf "USAGE:\n"
  printf "oncopmnet.sh <analysis-name> <input-bam-file> <nThreads> <input-hotspots-bed-file>\n\n"
}

exitWithError () {
  printf "\n***ERROR: $1\n\n"
  usage
  exit
}


if [[ ! -f "$INPUT_BAM" ]]; then
	exitWithError "Input BAM file $INPUT_BAM not found. Please check analysis name and/or bam file name."
fi

if [[ ! -f "$HOTSPOTS_BED" ]]; then
	exitWithError "Input Hotspots BED file $HOTSPOTS_BED not found."
fi

mkdir -p "$MOUNT_OUTDIR"

ALIGNED_BAM="$MOUNT_OUTDIR/${PREFIX}_aligned.bam"
SORTED_BAM="$MOUNT_OUTDIR/${PREFIX}_aligned_sorted.bam"

# The folowing commands (and especially their parameters) have been extracted
# from logs of runs in the Torrent Suite and documentation foujnd in the Thermo
# Fisher website and the web. They should be valid for DNA Ampliseq workflows,
# and more specifically the Oncomine Comprehensive Assay v3 and Colon and Lung
# designs.

BED_FNAME=$(basename "$TARGETS_BED" | sed 's/\.[^.]*$//')
PTRIM_BED="$MOUNT_OUTDIR/${BED_FNAME}_unmerged_detail.bed"
PLAIN_BED="$MOUNT_OUTDIR/${BED_FNAME}_merged_plain.bed"

printToLog "Preraring merged target BED file"
tvcutils validate_bed \
	--target-regions-bed "$TARGETS_BED" \
	--reference "$HGREF" \
	--unmerged-detail-bed "$PTRIM_BED" \
	--merged-plain-bed "$PLAIN_BED" || exit $?


HOTSPOTS_FNAME=$(basename "$HOTSPOTS_BED" | sed 's/\.[^.]*$//')
HOTSPOTS_LEFT_BED="$MOUNT_OUTDIR/${HOTSPOTS_FNAME}_left_aligned.bed"
HOTSPOTS_VCF="$MOUNT_OUTDIR/${HOTSPOTS_FNAME}.vcf"

printToLog "Preparing Hotspots VCF"
tvcutils prepare_hotspots \
	--input-bed "$HOTSPOTS_BED" \
	--reference "$HGREF" \
	--left-alignment on  --allow-block-substitutions on \
	--output-bed "$HOTSPOTS_LEFT_BED" \
	--output-vcf "$HOTSPOTS_VCF" \
	--unmerged-bed "$PTRIM_BED"

printToLog "Running tmap..."
tmap mapall -J 25 --end-repair 15 --do-repeat-clip --context \
	-u -v --prefix-exclude 5 -Y \
	-r "$INPUT_BAM" -f "$HGREF" -o 2 -n $THREADS \
	-i bam -s "$ALIGNED_BAM" stage1 map4 || exit $?

printToLog "Sorting aligned BAM..."
samtools sort -@ $THREADS -o "$SORTED_BAM" -O bam \
	-T `basename "$ALIGNED_BAM"` "$ALIGNED_BAM" || exit $?

printToLog "Indexing sorted BAM..."
samtools index -b "$SORTED_BAM" || exit $?

printToLog "Generating alignment stats..."
READS_RAW=`samtools stats "$INPUT_BAM" | grep '^SN' | grep 'raw total sequences:' | cut -f3`
READS_MAPPED=`samtools stats "$INPUT_BAM" | grep '^SN' | grep 'reads mapped:' | cut -f3`
READS_MAPQ20=`samtools view -c -q 20 "$INPUT_BAM"`
READ_LENGTH=`samtools stats "$INPUT_BAM" | grep '^SN' | grep 'average length:' | cut -f3`
echo -e "$READS_RAW \t $READS_MAPPED \t $READS_MAPQ20 \t $READ_LENGTH" >"$MOUNT_OUTDIR/${PREFIX}_alignment_stats.txt"

printToLog "Running variant caller..."
variant_caller_pipeline.py -o "$MOUNT_OUTDIR" \
    --num-threads $THREADS \
    --region-bed "$PLAIN_BED" \
    --primer-trim-bed "$PTRIM_BED" \
    --hotspot-vcf "$HOTSPOTS_VCF" \
    --input-bam "$SORTED_BAM" --parameters-file "$PARAMS_FILE" \
    --reference-fasta "$HGREF" || exit $?

chown -R 1004:1004 "$MOUNT_OUTDIR"