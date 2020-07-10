# Merged and Unmerged, Plain and Detail BED files

When "target" regions BED file is uploaded to torrent suite servers, it is 
converted into 4 different versions: merged and unmerged, plain and detailed.

* The merged plain version has the duplicate regions merged or removed. Use 
  the merged plain BED file with the TVC `--region-bed` option to avoid 
  redundant variant calls being made.

* Use the unmerged detail version of the BED file with the TVC 
  `--primer-trim-bed` option.

## How the merged/unmerged BED files are created

..using tvcutils validate_bed. The usage info reads:

```
Usage:   tvcutils validate_bed [options]

Input selection options (must provide one):
     --target-regions-bed        FILE       input is a target regions BED file (required)
     --hotspots-bed              FILE       input is a hotspots BED file (required)

General options:
     --reference                 FILE       FASTA file containing reference genome (required)
     --validation-log            FILE       log file for user-readable warning/error messages [stdout]
     --meta-json                 FILE       save validation and file statistics to json file [none]
     --unmerged-detail-bed       FILE       output a valid unmerged BED. To be used as input to --primer-trim-bed argument of variant_caller_pipeline.py (recommended) [none]
     --unmerged-plain-bed        FILE       output a valid unmerged BED. To be used as input to --region-bed argument of variant_caller_pipeline.py (recommended) [none]
     --merged-detail-bed         FILE       output an (almost) valid bedDetail merged BED [none]
     --merged-plain-bed          FILE       output a valid plain merged BED [none]
     --effective-bed             FILE       output a valid effective BED [none]
     --strict-check              on/off     exit 1 when there is error/filtered line [on]
```

So to produce all versions of an input targets ("designed") BED, we can use the following:

```bash
tvcutils validate_bed \
	--target-regions-bed Designed.bed \
	--reference hg19.fasta \
	--unmerged-detail-bed Designed.unmerged.detail.bed \
	--unmerged-plain-bed Designed.unmerged.plain.bed \
	--merged-detail-bed Designed.merged.detail.bed \
	--merged-plain-bed Designed.merged.plain.bed
```

It seems that TVC needs the "merged-plain" version as the target (region) BED..

