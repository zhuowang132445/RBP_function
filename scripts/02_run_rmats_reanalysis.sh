#!/bin/bash
set -euo pipefail

GROUP_DIR="${1:-as_reanalysis/grouping_balanced}"
OUT_BASE="${2:-as_reanalysis/rmats_recalc}"
RMATS_EXEC="${RMATS_EXEC:-$HOME/rmats-turbo-master/rmats.py}"
GTF_FILE="${GTF_FILE:-rice_all_genomes_v7.gtf}"
READ_LEN="${READ_LEN:-150}"
THREADS="${THREADS:-16}"
LIBTYPE="${LIBTYPE:-fr-unstranded}"
ALL_MANIFEST="${3:-as_reanalysis/metadata/bam_metrics_by_sample.tsv}"

mkdir -p "$OUT_BASE"

PROTEINS=(w1 w2 w3 w4 w5 w6)

for prot in "${PROTEINS[@]}"; do
  out_dir="${OUT_BASE}/${prot}"
  tmp_dir="${out_dir}/tmp"
  mkdir -p "$out_dir" "$tmp_dir"
  python "$RMATS_EXEC" \
    --b1 "${GROUP_DIR}/path_${prot}_high.txt" \
    --b2 "${GROUP_DIR}/path_${prot}_low.txt" \
    --gtf "$GTF_FILE" \
    -t paired \
    --readLength "$READ_LEN" \
    --nthread "$THREADS" \
    --od "$out_dir" \
    --tmp "$tmp_dir" \
    --libType "$LIBTYPE" \
    --variable-read-length \
    --allow-clipping
done

if [[ -f "$ALL_MANIFEST" ]]; then
  mkdir -p "${OUT_BASE}/all275_counts" "${OUT_BASE}/all275_counts/tmp"
  {
    printf "Sample\tBAM\n"
    awk -F'\t' 'NR>1 && $4 != "" {print $1 "\t" $4}' "$ALL_MANIFEST"
  } > "${OUT_BASE}/all275_counts/sample_bam_manifest.tsv"
  data_lines=$(( $(wc -l < "${OUT_BASE}/all275_counts/sample_bam_manifest.tsv") - 1 ))
  half=$(( data_lines / 2 ))
  tail -n +2 "${OUT_BASE}/all275_counts/sample_bam_manifest.tsv" | head -n "$half" | cut -f2 | paste -sd, - > "${OUT_BASE}/all275_counts/all_samples_part1.txt"
  tail -n +2 "${OUT_BASE}/all275_counts/sample_bam_manifest.tsv" | tail -n +"$((half + 1))" | cut -f2 | paste -sd, - > "${OUT_BASE}/all275_counts/all_samples_part2.txt"

  python "$RMATS_EXEC" \
    --b1 "${OUT_BASE}/all275_counts/all_samples_part1.txt" \
    --b2 "${OUT_BASE}/all275_counts/all_samples_part2.txt" \
    --gtf "$GTF_FILE" \
    -t paired \
    --readLength "$READ_LEN" \
    --nthread "$THREADS" \
    --od "${OUT_BASE}/all275_counts" \
    --tmp "${OUT_BASE}/all275_counts/tmp" \
    --libType "$LIBTYPE" \
    --variable-read-length \
    --allow-clipping
fi
