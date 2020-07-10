class: Workflow
cwlVersion: v1.0
id: torrent_pipeline_hotspots
label: torrent-pipeline-with-hotspots
$namespaces:
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: reference
    type: File
    doc: FASTA file containing reference genome
    secondaryFiles:
      - .fai
      - .tmap.anno
      - .tmap.bwt
      - .tmap.pac
      - .tmap.sa
    'sbg:x': -807.0403442382812
    'sbg:y': -457.4596862792969
  - id: parameters
    type: File
    'sbg:x': -810
    'sbg:y': -675
  - id: bed_file
    type: File
    'sbg:x': -801
    'sbg:y': 138
  - id: input_ubam
    type: File
    'sbg:x': -798
    'sbg:y': -99
  - id: hotspots_bed
    type: File
    'sbg:x': -797.9435424804688
    'sbg:y': -246.73387145996094
outputs:
  - id: coverage_report
    outputSource: coverage_analysis/coverage_report
    type: File
    'sbg:x': 415
    'sbg:y': 110
  - id: vcf
    outputSource: torrent_variant_caller/vcf
    type: File
    'sbg:x': 417
    'sbg:y': -706
  - id: small_variants_filtered
    outputSource: torrent_variant_caller/small_variants_filtered
    type: File
    'sbg:x': 412
    'sbg:y': -537
  - id: small_variants
    outputSource: torrent_variant_caller/small_variants
    type: File
    'sbg:x': 412
    'sbg:y': -359
  - id: indel_assembly
    outputSource: torrent_variant_caller/indel_assembly
    type: File
    'sbg:x': 419
    'sbg:y': -206
  - id: aligned_bam
    outputSource: samtools_index/alignments_with_index
    type: File
    'sbg:x': 411
    'sbg:y': -40
steps:
  - id: tmap
    in:
      - id: input_ubam
        source: input_ubam
      - id: reference
        source: reference
    out:
      - id: aligned_bam
    run: ../tools/tmap_ampliseq.cwl
    label: tmap
    'sbg:x': -445
    'sbg:y': -86
  - id: samtools_sort
    in:
      - id: input_bam
        source: tmap/aligned_bam
    out:
      - id: sorted_bam
    run: ../tools/samtools-sort.cwl
    'sbg:x': -299
    'sbg:y': -242
  - id: samtools_index
    in:
      - id: alignments
        source: samtools_sort/sorted_bam
    out:
      - id: alignments_with_index
    run: ../tools/samtools-index.cwl
    'sbg:x': -130
    'sbg:y': -283
  - id: coverage_analysis
    in:
      - id: reference
        source: reference
      - id: bed_file
        source: bed_file
      - id: aligned_bam
        source: samtools_index/alignments_with_index
    out:
      - id: coverage_report
    run: ../tools/cov-analysis.cwl
    'sbg:x': 86
    'sbg:y': 3
  - id: prepare_bed
    in:
      - id: reference
        source: reference
      - id: designed_bed
        source: bed_file
    out:
      - id: merged_plain_bed
      - id: unmerged_detail_bed
    run: ../tools/prepare_bed.cwl
    'sbg:x': -332.1209716796875
    'sbg:y': -507.01611328125
  - id: torrent_variant_caller
    in:
      - id: input_bam
        source: samtools_index/alignments_with_index
      - id: reference
        source: reference
      - id: plain_bed_file
        source: prepare_bed/merged_plain_bed
      - id: parameters
        source: parameters
      - id: ptrim_bed
        source: prepare_bed/unmerged_detail_bed
      - id: hotspots_vcf
        source: prepare_hotspots/hotspots_vcf
    out:
      - id: vcf
      - id: small_variants
      - id: indel_assembly
      - id: small_variants_filtered
    run: ../tools/torrent_variant_caller_hotspots.cwl
    label: torrent_variant_caller
    'sbg:x': 49.08064651489258
    'sbg:y': -649.1209716796875
  - id: prepare_hotspots
    in:
      - id: reference
        source: reference
      - id: hotspots_bed
        source: hotspots_bed
      - id: targets_unmerged_bed
        source: prepare_bed/unmerged_detail_bed
    out:
      - id: hotspots_vcf
      - id: hotspots_left_aligned_bed
    run: ../tools/prepare_hotspots.cwl
    'sbg:x': -168.21775817871094
    'sbg:y': -404.8064270019531
requirements: []
'sbg:toolAuthor': Stelios Sfakianakis
