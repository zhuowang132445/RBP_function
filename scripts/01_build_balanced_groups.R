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

expr_rdata <- get_arg("--expr_rdata", "salmon_cdna_mirna200_tpm_filtered_log_tpmc3_rankn.RData")
expr_object <- get_arg("--expr_object", "tpmc_mn_3")
gene_list_file <- get_arg("--gene_list", "RBP_gene_list.txt")
mapping_file <- get_arg("--mapping", "mapping_table.txt")
bam_metrics_file <- get_arg("--bam_metrics", "as_reanalysis/metadata/bam_metrics_by_sample.tsv")
manual_cov_file <- get_arg("--manual_covariates", "as_reanalysis/metadata/manual_covariates.tsv")
out_dir <- get_arg("--out_dir", "as_reanalysis/grouping_balanced")
group_frac <- as.numeric(get_arg("--group_frac", "0.25"))
latent_pcs <- as.integer(get_arg("--latent_pcs", "5"))
seed <- as.integer(get_arg("--seed", "20260423"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sort_fdbr <- function(x) {
  x <- unique(trimws(x))
  x <- x[x != ""]
  num <- suppressWarnings(as.integer(sub("^FDBR", "", x)))
  x[order(num, x, na.last = TRUE)]
}

read_mapping <- function(path) {
  df <- read.table(path, header = FALSE, skip = 1, sep = "", fill = TRUE, stringsAsFactors = FALSE)
  if (ncol(df) < 2) {
    stop("mapping_table.txt should contain at least 2 columns after the header.")
  }
  df <- df[, 1:2]
  colnames(df) <- c("FDBR", "Sample")
  df$FDBR <- trimws(df$FDBR)
  df$Sample <- trimws(df$Sample)
  unique(df[df$FDBR != "" & df$Sample != "", ])
}

compute_numeric_balance <- function(df, covars, group_col = "Group") {
  out <- list()
  if (length(covars) == 0) {
    return(data.frame())
  }
  for (covar in covars) {
    x1 <- df[df[[group_col]] == "high", covar]
    x0 <- df[df[[group_col]] == "low", covar]
    if (all(is.na(x1)) || all(is.na(x0))) {
      next
    }
    m1 <- mean(x1, na.rm = TRUE)
    m0 <- mean(x0, na.rm = TRUE)
    s1 <- sd(x1, na.rm = TRUE)
    s0 <- sd(x0, na.rm = TRUE)
    pooled <- sqrt((s1 ^ 2 + s0 ^ 2) / 2)
    smd <- ifelse(is.na(pooled) || pooled == 0, NA_real_, (m1 - m0) / pooled)
    pval <- tryCatch(t.test(x1, x0)$p.value, error = function(e) NA_real_)
    out[[covar]] <- data.frame(
      Covariate = covar,
      Type = "numeric",
      HighMean = m1,
      LowMean = m0,
      SMD = smd,
      PValue = pval
    )
  }
  do.call(rbind, out)
}

compute_factor_balance <- function(df, covars, group_col = "Group") {
  out <- list()
  if (length(covars) == 0) {
    return(data.frame())
  }
  for (covar in covars) {
    tab <- table(df[[group_col]], df[[covar]], useNA = "no")
    if (nrow(tab) < 2 || ncol(tab) < 1) {
      next
    }
    chi_p <- tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_)
    props <- prop.table(tab, 1)
    levels_cov <- colnames(tab)
    for (lev in levels_cov) {
      p1 <- props["high", lev]
      p0 <- props["low", lev]
      denom <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
      smd <- ifelse(is.na(denom) || denom == 0, NA_real_, (p1 - p0) / denom)
      out[[paste(covar, lev, sep = "::")]] <- data.frame(
        Covariate = covar,
        Type = paste0("factor:", lev),
        HighMean = p1,
        LowMean = p0,
        SMD = smd,
        PValue = chi_p
      )
    }
  }
  do.call(rbind, out)
}

choose_covariates <- function(df) {
  numeric_covars <- names(df)[vapply(df, is.numeric, logical(1))]
  factor_covars <- names(df)[vapply(df, function(x) is.character(x) || is.factor(x), logical(1))]
  factor_covars <- setdiff(factor_covars, c("Sample", "FDBR", "CRR", "BAM"))
  list(
    numeric = setdiff(numeric_covars, c("Expression", "propensity_score")),
    factor = setdiff(factor_covars, c("Group", "match_id"))
  )
}

greedy_match <- function(candidates, covariate_cols, exact_col = "SampleClass", seed = 1) {
  set.seed(seed)
  high <- candidates[candidates$Group == "high", , drop = FALSE]
  low <- candidates[candidates$Group == "low", , drop = FALSE]
  if (nrow(high) == 0 || nrow(low) == 0) {
    return(data.frame())
  }
  high <- high[order(high$propensity_score, decreasing = TRUE), ]
  used_low <- rep(FALSE, nrow(low))
  matches <- list()
  match_id <- 0
  caliper <- 0.5 * sd(qlogis(pmin(pmax(candidates$propensity_score, 1e-5), 1 - 1e-5)), na.rm = TRUE)
  if (is.na(caliper) || caliper <= 0) {
    caliper <- Inf
  }
  scaled <- scale(candidates[, covariate_cols, drop = FALSE])
  rownames(scaled) <- candidates$Sample
  for (i in seq_len(nrow(high))) {
    h <- high[i, ]
    pool_idx <- which(!used_low)
    if (!is.null(exact_col) && exact_col %in% names(low)) {
      pool_idx <- pool_idx[low[pool_idx, exact_col] == h[[exact_col]]]
    }
    if (length(pool_idx) == 0) {
      next
    }
    ps_diff <- abs(low$propensity_score[pool_idx] - h$propensity_score)
    pool_idx <- pool_idx[ps_diff <= caliper]
    if (length(pool_idx) == 0) {
      next
    }
    hvec <- scaled[h$Sample, , drop = FALSE]
    lmat <- scaled[low$Sample[pool_idx], , drop = FALSE]
    dists <- rowSums((t(t(lmat) - as.numeric(hvec))) ^ 2, na.rm = TRUE)
    best <- pool_idx[which.min(dists)]
    match_id <- match_id + 1
    h$match_id <- match_id
    low_row <- low[best, , drop = FALSE]
    low_row$match_id <- match_id
    matches[[length(matches) + 1]] <- h
    matches[[length(matches) + 1]] <- low_row
    used_low[best] <- TRUE
  }
  if (length(matches) == 0) {
    return(data.frame())
  }
  do.call(rbind, matches)
}

load(expr_rdata)
if (!exists(expr_object)) {
  stop(sprintf("Object %s not found in %s", expr_object, expr_rdata))
}
expr_mat <- get(expr_object)
colnames(expr_mat) <- trimws(colnames(expr_mat))
rownames(expr_mat) <- trimws(rownames(expr_mat))

mapping_df <- read_mapping(mapping_file)
sample_to_fdbr <- setNames(mapping_df$FDBR, mapping_df$Sample)

metrics_df <- NULL
if (file.exists(bam_metrics_file)) {
  metrics_df <- read.delim(bam_metrics_file, header = TRUE, check.names = FALSE)
}

if (is.null(metrics_df) || nrow(metrics_df) == 0) {
  fdbr_ordered <- sort_fdbr(mapping_df$FDBR)
  bam_files <- list.files("/public/home/zhuowang/smoke/egwas/02_StarAlign",
                          pattern = "_new\\.Aligned\\.sortedByCoord\\.out\\.bam$",
                          full.names = TRUE)
  bam_files <- sort(bam_files)
  crr_ids <- sub("_new\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(bam_files))
  metrics_df <- data.frame(
    Sample = mapping_df$Sample[match(fdbr_ordered, mapping_df$FDBR)],
    FDBR = fdbr_ordered,
    CRR = crr_ids,
    BAM = bam_files,
    mapped_reads = NA_real_,
    total_reads = NA_real_,
    batch = "unknown",
    stringsAsFactors = FALSE
  )
}

metrics_df$Sample <- trimws(metrics_df$Sample)
metrics_df$FDBR <- trimws(metrics_df$FDBR)
metrics_df$SampleClass <- substr(metrics_df$Sample, 1, 1)
metrics_df$log_mapped_reads <- log10(as.numeric(metrics_df$mapped_reads) + 1)
metrics_df$batch <- ifelse(is.na(metrics_df$batch) | metrics_df$batch == "", "unknown", metrics_df$batch)

if (file.exists(manual_cov_file)) {
  manual_df <- read.delim(manual_cov_file, header = TRUE, check.names = FALSE)
  manual_df$Sample <- trimws(manual_df$Sample)
  metrics_df <- merge(metrics_df, manual_df, by = "Sample", all.x = TRUE, sort = FALSE)
}

gene_df <- read.table(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
gene_df <- gene_df[, 1:2]
colnames(gene_df) <- c("GeneID", "Protein")
gene_df$GeneID <- trimws(gene_df$GeneID)
gene_df$Protein <- trimws(gene_df$Protein)

all_rbp_ids <- unique(gene_df$GeneID)
expr_for_pca <- expr_mat[setdiff(rownames(expr_mat), all_rbp_ids), , drop = FALSE]
expr_for_pca <- expr_for_pca[, intersect(colnames(expr_for_pca), metrics_df$Sample), drop = FALSE]
expr_for_pca <- expr_for_pca[apply(expr_for_pca, 1, function(x) sd(as.numeric(x), na.rm = TRUE) > 0), , drop = FALSE]
expr_for_pca <- t(scale(t(expr_for_pca)))
expr_for_pca[is.na(expr_for_pca)] <- 0
pca <- prcomp(t(expr_for_pca), center = FALSE, scale. = FALSE)
pc_keep <- seq_len(min(latent_pcs, ncol(pca$x)))
pc_df <- data.frame(Sample = rownames(pca$x), pca$x[, pc_keep, drop = FALSE], check.names = FALSE)
colnames(pc_df)[-1] <- paste0("ExprPC", pc_keep)

cov_df <- merge(metrics_df, pc_df, by = "Sample", all.x = TRUE, sort = FALSE)

for (nm in names(cov_df)) {
  if (is.numeric(cov_df[[nm]]) && anyNA(cov_df[[nm]])) {
    med <- median(cov_df[[nm]], na.rm = TRUE)
    if (is.finite(med)) {
      cov_df[[nm]][is.na(cov_df[[nm]])] <- med
    }
  }
  if ((is.character(cov_df[[nm]]) || is.factor(cov_df[[nm]])) && anyNA(cov_df[[nm]])) {
    cov_df[[nm]][is.na(cov_df[[nm]])] <- "unknown"
  }
}

write.table(cov_df, file.path(out_dir, "sample_covariates_used.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_lines <- c()

for (i in seq_len(nrow(gene_df))) {
  gene_id <- gene_df$GeneID[i]
  prot <- gene_df$Protein[i]
  if (!(gene_id %in% rownames(expr_mat))) {
    next
  }

  expr_vec <- as.numeric(expr_mat[gene_id, ])
  names(expr_vec) <- colnames(expr_mat)
  sample_df <- cov_df[cov_df$Sample %in% names(expr_vec), , drop = FALSE]
  sample_df$Expression <- expr_vec[sample_df$Sample]
  sample_df <- sample_df[!is.na(sample_df$Expression), , drop = FALSE]
  sample_df <- sample_df[order(sample_df$Expression, decreasing = TRUE), ]
  n_group <- floor(nrow(sample_df) * group_frac)
  high_candidates <- sample_df[seq_len(n_group), , drop = FALSE]
  low_candidates <- sample_df[(nrow(sample_df) - n_group + 1):nrow(sample_df), , drop = FALSE]
  high_candidates$Group <- "high"
  low_candidates$Group <- "low"
  candidates <- rbind(high_candidates, low_candidates)
  candidates$batch <- ifelse(is.na(candidates$batch) | candidates$batch == "", "unknown", candidates$batch)
  candidates$SampleClass <- ifelse(is.na(candidates$SampleClass), "U", candidates$SampleClass)

  cov_choices <- choose_covariates(candidates)
  num_covars <- cov_choices$numeric
  fac_covars <- cov_choices$factor
  num_covars <- setdiff(num_covars, c("mapped_reads", "total_reads"))
  if ("log_mapped_reads" %in% names(candidates)) {
    num_covars <- union(num_covars, "log_mapped_reads")
  }
  model_terms <- c(num_covars, fac_covars)
  model_terms <- model_terms[model_terms %in% names(candidates)]
  candidates$GroupBinary <- ifelse(candidates$Group == "high", 1, 0)
  if (length(model_terms) == 0) {
    candidates$propensity_score <- 0.5
  } else {
    formula_str <- paste("GroupBinary ~", paste(model_terms, collapse = " + "))
    fit <- tryCatch(glm(as.formula(formula_str), data = candidates, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) {
      candidates$propensity_score <- 0.5
    } else {
      candidates$propensity_score <- pmin(pmax(fitted(fit), 1e-5), 1 - 1e-5)
    }
  }

  match_covars <- num_covars[num_covars %in% names(candidates)]
  if (length(match_covars) == 0) {
    match_covars <- "Expression"
  }
  matched <- greedy_match(candidates, covariate_cols = match_covars, exact_col = "SampleClass", seed = seed)
  if (nrow(matched) == 0) {
    warning(sprintf("No balanced matches found for %s", prot))
    next
  }

  matched <- matched[order(matched$match_id, ifelse(matched$Group == "high", 0, 1)), ]
  matched$RBP <- prot
  matched$GeneID <- gene_id

  pre_num <- compute_numeric_balance(candidates, num_covars)
  pre_fac <- compute_factor_balance(candidates, fac_covars)
  post_num <- compute_numeric_balance(matched, num_covars)
  post_fac <- compute_factor_balance(matched, fac_covars)
  pre_bal <- rbind(pre_num, pre_fac)
  post_bal <- rbind(post_num, post_fac)
  write.table(pre_bal, file.path(out_dir, paste0(prot, "_balance_pre.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(post_bal, file.path(out_dir, paste0(prot, "_balance_post.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(matched, file.path(out_dir, paste0(prot, "_grouping_balanced.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  high_paths <- matched$BAM[matched$Group == "high"]
  low_paths <- matched$BAM[matched$Group == "low"]
  writeLines(paste(high_paths, collapse = ","), file.path(out_dir, paste0("path_", prot, "_high.txt")))
  writeLines(paste(low_paths, collapse = ","), file.path(out_dir, paste0("path_", prot, "_low.txt")))
  writeLines(matched$Sample[matched$Group == "high"], file.path(out_dir, paste0(prot, "_high_samples.txt")))
  writeLines(matched$Sample[matched$Group == "low"], file.path(out_dir, paste0(prot, "_low_samples.txt")))

  summary_lines <- c(
    summary_lines,
    sprintf("%s\tgene=%s\thigh=%d\tlow=%d\tcandidate_high=%d\tcandidate_low=%d",
            prot, gene_id, sum(matched$Group == "high"), sum(matched$Group == "low"),
            sum(candidates$Group == "high"), sum(candidates$Group == "low"))
  )
}

writeLines(summary_lines, file.path(out_dir, "grouping_summary.txt"))
