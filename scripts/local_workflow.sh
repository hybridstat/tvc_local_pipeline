#!/bin/bash

printToLog() {
  date +"%F__%T ---> $1"
}

usage () {
  printf "USAGE:\n"
  printf "oncopmnet.sh <analysis-name> <input-bam-file> <nThreads>\n\n"
}

exitWithError () {
  printf "\n***ERROR: $1\n\n"
  usage
  exit 1
}

if [[ ! -f "$INPUT_BAM" ]]; then
  exitWithError "Input BAM file $INPUT_BAM not found. Please check analysis name and/or bam file name."
fi

# if [ -d "$MOUNT_OUTDIR" ]; then
#     exitWithError "Directory $MOUNT_OUTDIR exists. Sample $PREFIX has already been analyzed under analysis $ANALYSIS_NAME. Please remove $PREFIX subdirectory or choose a different analysis name or BAM sample for this analysis."
# fi

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

printToLog "Preraring merged (plain) target BED file..."
tvcutils validate_bed \
	--target-regions-bed "$TARGETS_BED" \
	--reference "$HGREF" \
	--unmerged-detail-bed "$PTRIM_BED" \
	--merged-plain-bed "$PLAIN_BED" || exit $?

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
    --input-bam "$SORTED_BAM" --parameters-file "$PARAMS_FILE" \
    --reference-fasta "$HGREF" || exit $?

chown -R 1004:1004 "$MOUNT_OUTDIR"