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

rmats_base <- get_arg("--rmats_base", "as_reanalysis/rmats_recalc")
cont_dir <- get_arg("--cont_dir", "as_reanalysis/continuous_model")
out_dir <- get_arg("--out_dir", "as_reanalysis/integration")
hl_fdr <- as.numeric(get_arg("--hl_fdr", "0.05"))
hl_delta <- as.numeric(get_arg("--hl_delta", "0.10"))
cont_fdr <- as.numeric(get_arg("--cont_fdr", "0.05"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

proteins <- c("w1", "w2", "w3", "w4", "w5", "w6")
event_types <- c("SE", "RI", "A3SS", "A5SS", "MXE")

hl_all <- list()
for (prot in proteins) {
  for (ev in event_types) {
    path <- file.path(rmats_base, prot, paste0(ev, ".MATS.JCEC.txt"))
    if (!file.exists(path)) {
      next
    }
    df <- read.delim(path, header = TRUE, check.names = FALSE)
    df$EventUID <- paste(ev, df$ID, sep = "|")
    df$RBP <- prot
    keep <- c("RBP", "EventUID", "GeneID", "geneSymbol", "FDR", "IncLevelDifference")
    keep <- keep[keep %in% names(df)]
    out <- df[, keep, drop = FALSE]
    names(out)[names(out) == "FDR"] <- "HL_FDR"
    names(out)[names(out) == "IncLevelDifference"] <- "HL_DeltaPSI"
    hl_all[[paste(prot, ev, sep = "_")]] <- out
  }
}

hl_df <- if (length(hl_all) > 0) do.call(rbind, hl_all) else data.frame()
cont_files <- list.files(cont_dir, pattern = "_continuous_model\\.tsv$", full.names = TRUE)
cont_df <- if (length(cont_files) > 0) {
  do.call(rbind, lapply(cont_files, function(path) read.delim(path, header = TRUE, check.names = FALSE)))
} else {
  data.frame()
}

if (nrow(hl_df) == 0 && nrow(cont_df) == 0) {
  merged <- data.frame()
} else if (nrow(hl_df) == 0) {
  merged <- cont_df
} else if (nrow(cont_df) == 0) {
  merged <- hl_df
} else {
  merged <- merge(hl_df, cont_df, by = c("RBP", "EventUID"), all = TRUE, sort = FALSE)
}

pick_first_nonempty <- function(a, b) {
  a <- trimws(as.character(a))
  b <- trimws(as.character(b))
  a[is.na(a)] <- ""
  b[is.na(b)] <- ""
  out <- ifelse(a != "", a, b)
  out[out == ""] <- NA_character_
  out
}

if ("GeneID.x" %in% names(merged) && "GeneID.y" %in% names(merged)) {
  merged$GeneID <- pick_first_nonempty(merged$GeneID.x, merged$GeneID.y)
  merged$GeneID.x <- NULL
  merged$GeneID.y <- NULL
}
if ("geneSymbol.x" %in% names(merged) && "geneSymbol.y" %in% names(merged)) {
  merged$geneSymbol <- pick_first_nonempty(merged$geneSymbol.x, merged$geneSymbol.y)
  merged$geneSymbol.x <- NULL
  merged$geneSymbol.y <- NULL
}

if (!("HL_FDR" %in% names(merged))) merged$HL_FDR <- NA_real_
if (!("HL_DeltaPSI" %in% names(merged))) merged$HL_DeltaPSI <- NA_real_
if (!("FDR" %in% names(merged))) merged$FDR <- NA_real_
if (!("Beta" %in% names(merged))) merged$Beta <- NA_real_
merged$HL_Significant <- !is.na(merged$HL_FDR) & merged$HL_FDR < hl_fdr & abs(merged$HL_DeltaPSI) >= hl_delta
merged$Continuous_Significant <- !is.na(merged$FDR) & merged$FDR < cont_fdr
merged$Direction_Consistent <- !is.na(merged$HL_DeltaPSI) & !is.na(merged$Beta) &
  sign(merged$HL_DeltaPSI) == sign(merged$Beta)
merged$Consensus <- merged$HL_Significant & merged$Continuous_Significant & merged$Direction_Consistent

write.table(merged, file.path(out_dir, "consensus_events.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

consensus_only <- merged[merged$Consensus, , drop = FALSE]
if (nrow(consensus_only) > 0) {
  genes <- aggregate(EventUID ~ RBP + GeneID + geneSymbol, data = consensus_only, FUN = length)
  names(genes)[names(genes) == "EventUID"] <- "ConsensusEventCount"
  min_hl <- aggregate(HL_FDR ~ RBP + GeneID + geneSymbol, data = consensus_only, FUN = min)
  min_cont <- aggregate(FDR ~ RBP + GeneID + geneSymbol, data = consensus_only, FUN = min)
  genes <- Reduce(function(x, y) merge(x, y, by = c("RBP", "GeneID", "geneSymbol"), all = TRUE),
                  list(genes, min_hl, min_cont))
} else {
  genes <- data.frame()
}

write.table(genes, file.path(out_dir, "consensus_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (nrow(merged) > 0) {
  summary_df <- aggregate(
    cbind(HL_Significant, Continuous_Significant, Direction_Consistent, Consensus) ~ RBP,
    data = transform(merged,
      HL_Significant = as.integer(HL_Significant),
      Continuous_Significant = as.integer(Continuous_Significant),
      Direction_Consistent = as.integer(Direction_Consistent),
      Consensus = as.integer(Consensus)
    ),
    FUN = sum
  )
  write.table(summary_df, file.path(out_dir, "summary_by_rbp.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
