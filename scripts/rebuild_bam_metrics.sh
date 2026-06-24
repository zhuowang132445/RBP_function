#!/bin/bash
set -euo pipefail

mkdir -p as_reanalysis/metadata

BAM_DIR="/public/home/zhuowang/smoke/egwas/02_StarAlign"

paste rmats_group_files/fdbr_ordered_from_mapping.txt rmats_group_files/bam_sorted_names.txt \
> as_reanalysis/metadata/fdbr_bam_basename.tsv

awk -v bam_dir="$BAM_DIR" 'BEGIN{FS=OFS="\t"}
NR==FNR{
  if(FNR==1) next
  fdbr=$1
  sample=$2
  gsub(/^[ \t]+|[ \t]+$/, "", fdbr)
  gsub(/^[ \t]+|[ \t]+$/, "", sample)
  if(fdbr!="" && sample!=""){
    sample_for_fdbr[fdbr]=sample
  }
  next
}
{
  fdbr=$1
  bamfile=$2
  crr=bamfile
  sub(/_new.Aligned.sortedByCoord.out.bam$/, "", crr)
  sample=sample_for_fdbr[fdbr]
  if(sample!=""){
    print sample, fdbr, crr, bam_dir "/" bamfile
  }
}' rmats_group_files/clean_sample_mapping.txt as_reanalysis/metadata/fdbr_bam_basename.tsv \
> as_reanalysis/metadata/sample_bam_base.tsv

{
  printf 'Sample\tFDBR\tCRR\tBAM\tmapped_reads\ttotal_reads\tbatch\n'
  while IFS=$'\t' read -r sample fdbr crr bam; do
    mapped=$(samtools idxstats "$bam" | awk '{m+=$3} END{print m+0}')
    total=$(samtools idxstats "$bam" | awk '{m+=$3+$4} END{print m+0}')
    batch=$(samtools view -H "$bam" | awk '
      BEGIN{batch="unknown"}
      /^@RG/{
        for(i=1;i<=NF;i++){
          if($i ~ /^PU:/){sub(/^PU:/,"",$i); batch=$i}
          else if($i ~ /^LB:/ && batch=="unknown"){sub(/^LB:/,"",$i); batch=$i}
          else if($i ~ /^ID:/ && batch=="unknown"){sub(/^ID:/,"",$i); batch=$i}
        }
      }
      END{print batch}
    ')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$fdbr" "$crr" "$bam" "$mapped" "$total" "$batch"
  done < as_reanalysis/metadata/sample_bam_base.tsv
} > as_reanalysis/metadata/bam_metrics_by_sample.tsv
