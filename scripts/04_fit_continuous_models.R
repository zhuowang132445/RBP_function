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
covariates_file <- get_arg("--covariates", "as_reanalysis/grouping_balanced/sample_covariates_used.tsv")
psi_file <- get_arg("--psi", "as_reanalysis/psi_matrix/events_psi_matrix.tsv")
out_dir <- get_arg("--out_dir", "as_reanalysis/continuous_model")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

clamp01 <- function(x) pmin(pmax(x, 1e-4), 1 - 1e-4)

load(expr_rdata)
expr_mat <- get(expr_object)
colnames(expr_mat) <- trimws(colnames(expr_mat))
rownames(expr_mat) <- trimws(rownames(expr_mat))

gene_df <- read.table(gene_list_file, header = FALSE, stringsAsFactors = FALSE)
gene_df <- gene_df[, 1:2]
colnames(gene_df) <- c("GeneID", "Protein")

cov_df <- read.delim(covariates_file, header = TRUE, check.names = FALSE)
psi_df <- read.delim(psi_file, header = TRUE, check.names = FALSE)

meta_cols <- c("EventUID", "EventType", "ID", "GeneID", "geneSymbol")
sample_cols <- setdiff(colnames(psi_df), meta_cols)

numeric_covars <- names(cov_df)[vapply(cov_df, is.numeric, logical(1))]
numeric_covars <- setdiff(numeric_covars, c("mapped_reads", "total_reads"))
if ("log_mapped_reads" %in% names(cov_df)) {
  numeric_covars <- union(numeric_covars, "log_mapped_reads")
}
numeric_covars <- setdiff(numeric_covars, "Expression")
factor_covars <- c()
if ("SampleClass" %in% names(cov_df) && length(unique(cov_df$SampleClass)) > 1) {
  factor_covars <- c(factor_covars, "SampleClass")
}
if ("batch" %in% names(cov_df) && length(unique(cov_df$batch)) > 1) {
  factor_covars <- c(factor_covars, "batch")
}

for (i in seq_len(nrow(gene_df))) {
  gene_id <- trimws(gene_df$GeneID[i])
  prot <- trimws(gene_df$Protein[i])
  if (!(gene_id %in% rownames(expr_mat))) {
    next
  }
  rbp_expr <- as.numeric(expr_mat[gene_id, ])
  names(rbp_expr) <- colnames(expr_mat)
  dat <- cov_df[cov_df$Sample %in% sample_cols, , drop = FALSE]
  dat$RBPExpression <- scale(rbp_expr[dat$Sample])[, 1]

  rhs <- c("RBPExpression", numeric_covars, factor_covars)
  rhs <- rhs[rhs %in% names(dat)]
  rhs <- unique(rhs)
  form <- as.formula(paste("PSI_logit ~", paste(rhs, collapse = " + ")))

  res_list <- vector("list", nrow(psi_df))
  for (j in seq_len(nrow(psi_df))) {
    event_vals <- as.numeric(psi_df[j, dat$Sample, drop = TRUE])
    work <- dat
    work$PSI <- event_vals
    work <- work[!is.na(work$PSI), , drop = FALSE]
    if (nrow(work) < 50) {
      next
    }
    if (sd(work$PSI, na.rm = TRUE) < 0.02) {
      next
    }
    work$PSI_logit <- qlogis(clamp01(work$PSI))
    fit <- tryCatch(lm(form, data = work), error = function(e) NULL)
    if (is.null(fit)) {
      next
    }
    coefs <- summary(fit)$coefficients
    if (!("RBPExpression" %in% rownames(coefs))) {
      next
    }
    res_list[[j]] <- data.frame(
      RBP = prot,
      RBP_GeneID = gene_id,
      EventUID = psi_df$EventUID[j],
      EventType = psi_df$EventType[j],
      GeneID = psi_df$GeneID[j],
      geneSymbol = psi_df$geneSymbol[j],
      N = nrow(work),
      Beta = coefs["RBPExpression", "Estimate"],
      SE = coefs["RBPExpression", "Std. Error"],
      T = coefs["RBPExpression", "t value"],
      PValue = coefs["RBPExpression", "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }
  res_df <- do.call(rbind, res_list)
  if (is.null(res_df) || nrow(res_df) == 0) {
    next
  }
  res_df$FDR <- p.adjust(res_df$PValue, method = "BH")
  write.table(res_df, file.path(out_dir, paste0(prot, "_continuous_model.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
