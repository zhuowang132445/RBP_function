suppressWarnings(suppressMessages({
  options(stringsAsFactors = FALSE)
}))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[idx + 1]
}

rmats_dir <- get_arg("--rmats_dir", "as_reanalysis/rmats_recalc/all275_counts")
manifest_file <- get_arg("--manifest", "as_reanalysis/rmats_recalc/all275_counts/sample_bam_manifest.tsv")
out_dir <- get_arg("--out_dir", "as_reanalysis/psi_matrix")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

detect_cols <- function(df) {
  inc1 <- intersect(c("IC_SAMPLE_1", "IJC_SAMPLE_1"), colnames(df))
  skip1 <- intersect(c("SC_SAMPLE_1", "SJC_SAMPLE_1"), colnames(df))
  inc2 <- intersect(c("IC_SAMPLE_2", "IJC_SAMPLE_2"), colnames(df))
  skip2 <- intersect(c("SC_SAMPLE_2", "SJC_SAMPLE_2"), colnames(df))
  if (length(inc1) == 0 || length(skip1) == 0 || length(inc2) == 0 || length(skip2) == 0) {
    stop("Could not detect count columns in rMATS output.")
  }
  list(inc1 = inc1[1], skip1 = skip1[1], inc2 = inc2[1], skip2 = skip2[1])
}

calc_psi <- function(inc_str, skip_str, inc_len, skip_len) {
  inc <- as.numeric(strsplit(inc_str, ",", fixed = TRUE)[[1]])
  skip <- as.numeric(strsplit(skip_str, ",", fixed = TRUE)[[1]])
  val <- (inc / inc_len) / ((inc / inc_len) + (skip / skip_len))
  val[is.nan(val) | is.infinite(val)] <- NA_real_
  val
}

sample_manifest <- read.delim(manifest_file, header = TRUE, check.names = FALSE)
all_samples <- sample_manifest$Sample
half <- floor(nrow(sample_manifest) / 2)
group1_samples <- sample_manifest$Sample[seq_len(half)]
group2_samples <- sample_manifest$Sample[(half + 1):nrow(sample_manifest)]

event_types <- c("SE", "RI", "A3SS", "A5SS", "MXE")
psi_list <- list()
meta_list <- list()

for (ev in event_types) {
  path <- file.path(rmats_dir, paste0(ev, ".MATS.JCEC.txt"))
  if (!file.exists(path)) {
    next
  }
  df <- read.delim(path, header = TRUE, check.names = FALSE)
  cols <- detect_cols(df)
  event_uid <- paste(ev, df$ID, sep = "|")
  psi_mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(all_samples),
                    dimnames = list(event_uid, all_samples))

  for (i in seq_len(nrow(df))) {
    psi1 <- calc_psi(df[[cols$inc1]][i], df[[cols$skip1]][i], df$IncFormLen[i], df$SkipFormLen[i])
    psi2 <- calc_psi(df[[cols$inc2]][i], df[[cols$skip2]][i], df$IncFormLen[i], df$SkipFormLen[i])
    psi_mat[i, group1_samples] <- psi1
    psi_mat[i, group2_samples] <- psi2
  }

  psi_list[[ev]] <- psi_mat
  meta_list[[ev]] <- data.frame(
    EventUID = event_uid,
    EventType = ev,
    ID = df$ID,
    GeneID = if ("GeneID" %in% names(df)) df$GeneID else NA,
    geneSymbol = if ("geneSymbol" %in% names(df)) df$geneSymbol else NA,
    stringsAsFactors = FALSE
  )
}

psi_all <- do.call(rbind, psi_list)
meta_all <- do.call(rbind, meta_list)

psi_out <- cbind(meta_all, as.data.frame(psi_all, check.names = FALSE))
write.table(psi_out, file.path(out_dir, "events_psi_matrix.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
