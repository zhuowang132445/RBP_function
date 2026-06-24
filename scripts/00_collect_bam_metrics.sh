#!/bin/bash
set -euo pipefail

OUT_DIR="${1:-as_reanalysis/metadata}"
MAPPING_FILE="${2:-mapping_table.txt}"
BAM_DIR="${3:-/public/home/zhuowang/smoke/egwas/02_StarAlign}"
BAM_SUFFIX="${4:-_new.Aligned.sortedByCoord.out.bam}"
THREADS="${THREADS:-4}"

mkdir -p "$OUT_DIR"

cat > "${OUT_DIR}/collect_bam_metrics.awk" <<'AWK'
BEGIN {
  OFS = "\t";
}
FNR == NR && NR > 1 {
  fdbr = $1;
  sample = $2;
  gsub(/^[ \t]+|[ \t]+$/, "", fdbr);
  gsub(/^[ \t]+|[ \t]+$/, "", sample);
  if (fdbr != "" && sample != "") {
    sample_for_fdbr[fdbr] = sample;
  }
  next;
}
FNR == 1 {
  next;
}
{
  fdbr = $1;
  crr = $2;
  bam = $3;
  mapped = $4;
  total = $5;
  batch = $6;
  sample = sample_for_fdbr[fdbr];
  if (sample == "") {
    sample = "NA";
  }
  print sample, fdbr, crr, bam, mapped, total, batch;
}
AWK

printf "Sample\tFDBR\n" > "${OUT_DIR}/mapping_table_clean.tsv"
tail -n +2 "$MAPPING_FILE" | awk 'NF >= 2 {print $1 "\t" $2}' >> "${OUT_DIR}/mapping_table_clean.tsv"

find "$BAM_DIR" -maxdepth 1 -name "*${BAM_SUFFIX}" | sort > "${OUT_DIR}/bam_paths_sorted.txt"

awk -v suffix="$BAM_SUFFIX" '
{
  bam = $0;
  file = $0;
  sub(/^.*\//, "", file);
  crr = file;
  sub(suffix "$", "", crr);
  print crr "\t" bam;
}
' "${OUT_DIR}/bam_paths_sorted.txt" > "${OUT_DIR}/crr_to_bam.tsv"

awk '{
  n++;
  fdbr = "FDBR" n;
  print fdbr "\t" $1 "\t" $2;
}' "${OUT_DIR}/crr_to_bam.tsv" > "${OUT_DIR}/fdbr_crr_bam.tsv"

{
  printf "FDBR\tCRR\tBAM\tmapped_reads\ttotal_reads\tbatch\n";
  while IFS=$'\t' read -r fdbr crr bam; do
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
    printf "%s\t%s\t%s\t%s\t%s\t%s\n", fdbr, crr, bam, mapped, total, batch;
  done < "${OUT_DIR}/fdbr_crr_bam.tsv"
} > "${OUT_DIR}/bam_metrics_fdbr.tsv"

awk -f "${OUT_DIR}/collect_bam_metrics.awk" \
  "${OUT_DIR}/mapping_table_clean.tsv" \
  "${OUT_DIR}/bam_metrics_fdbr.tsv" > "${OUT_DIR}/bam_metrics_by_sample.tsv"

rm -f "${OUT_DIR}/collect_bam_metrics.awk"

echo "Wrote:"
echo "  ${OUT_DIR}/bam_metrics_fdbr.tsv"
echo "  ${OUT_DIR}/bam_metrics_by_sample.tsv"
