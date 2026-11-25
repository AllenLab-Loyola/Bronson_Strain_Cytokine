## =====================================================================
## Figure S2: laboratory and clinical metrics by source-patient outcome
## Bronson et al., strain-dependent Pa bacteremia manuscript
##
##   1. Classifies each isolate by whether its source patient SURVIVED or
##      DIED during the index admission.
##   2. Compares each metric between the groups: day-0 SOFA, murine
##      virulence score, total cytokine induction, IL-1beta, IL-6, TNFalpha.
##   3. Draws a box-and-dot plot per comparison, combined into one figure
##      and also saved individually.
##
## USAGE
##   Rscript FigS2_mortality_comparison.R
##   Rscript FigS2_mortality_comparison.R --mortality X.xlsx --cytokine Y.xlsx
##
## INPUTS   patient.xlsx, cytokine.xlsx -- read with the same loader and
##          rules as Fig4_clinical_power_analysis.R, so both figures use
##          the identical 85-isolate cohort
## OUTPUTS  classified_isolates.csv     per-isolate group + values
##          mortality_comparison.csv    the statistics below
##          FigS2_mortality_comparison.png / .pdf   combined 2 x 3 figure
##          panels/<metric>.png         each comparison on its own
##
## DEPENDENCIES  readxl, ggplot2, patchwork; ragg recommended for Greek
##               glyphs in the PNGs (cairo is used otherwise)
##
## STATISTICS
##   Mann-Whitney is the primary test: Shapiro-Wilk rejects normality for
##   most of these measures, and the virulence score is an integer 1-7.
##   A Welch t-test is reported alongside for reference only.
##   Effect size is the rank-biserial correlation with a percentile
##   bootstrap CI; POSITIVE means higher values among non-survivors.
##   Holm adjustment covers the FIVE laboratory metrics. Day-0 SOFA is a
##   prespecified positive control, not one of the hypotheses under test,
##   so it is excluded from the correction (HOLM_INCLUDES_SOFA = TRUE
##   folds it in).
##
## INTERPRETING A NULL RESULT
##   These comparisons are underpowered for anything but a large effect:
##   with 31 deaths among 85 cases, 80% power requires |r| = 0.36 and
##   every laboratory metric here is at or below 0.17. A non-significant
##   p-value bounds the effect, it does not exclude one. See
##   Fig4_clinical_power_analysis.R for the power analysis.
## =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(patchwork)
})


## ---- Settings ---------------------------------------------------------

set.seed(20261119)
NBOOT <- 10000L
HOLM_INCLUDES_SOFA <- FALSE

## column, axis label, is-positive-control
VARS <- list(
  list(col = "SOFA",      lab = "Day-0 SOFA score",                  ctrl = TRUE),
  list(col = "Virulence", lab = "Murine virulence score",            ctrl = FALSE),
  list(col = "TOTAL",     lab = "Total cytokine\n(sum log10 pg/mL)", ctrl = FALSE),
  list(col = "IL_1B",     lab = "IL-1β (log10 pg/mL)",          ctrl = FALSE),
  list(col = "IL6",       lab = "IL-6 (log10 pg/mL)",                ctrl = FALSE),
  list(col = "TNF",       lab = "TNFα (log10 pg/mL)",           ctrl = FALSE))

COL_SURV <- "#2c6fbb"
COL_DIED <- "#c0392b"

## ---- Arguments --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
f_mort <- getarg("--mortality", "patient.xlsx")
f_cyto <- getarg("--cytokine",  "cytokine.xlsx")


## ---- Data -------------------------------------------------------------
## Same loader as Fig4_clinical_power_analysis.R: change one, change both.

## Canonical column names and the spelling variants seen in practice.
## Matching ignores case, spaces, underscores, hyphens and dots.
ALIASES <- list(
  Isolate   = c("isolate", "strain", "strainid", "isolateid"),
  Mortality = c("mortality", "died", "death", "dead",
                "inhospitalmortality", "mortality01"),
  SOFA      = c("sofa", "day0sofa", "sofaday0", "d0sofa", "sofad0",
                "sofascore", "day0sofascore"),
  IL_1B     = c("il1b", "il1beta", "il1b_pgml", "illb", "il1"),
  IL6       = c("il6", "il6_pgml"),
  TNF       = c("tnf", "tnfa", "tnfalpha", "tnfa_pgml", "tnf_pgml"),
  Virulence = c("virulence", "virulencescore", "virscore",
                "invivovirulence", "virulencegroup",
                "mousevirulence", "mousevirulencescore",
                "murinevirulence", "murinevirulencescore"),
  Replicate = c("replicate", "rep")
)

canon_key <- function(x) gsub("[^a-z0-9]", "", tolower(x))

canonicalize <- function(df, label) {
  orig <- names(df)
  keys <- canon_key(orig)
  renamed <- character(0)
  for (target in names(ALIASES)) {
    if (target %in% names(df)) next
    hit <- which(keys %in% ALIASES[[target]])
    if (length(hit) == 1L) {
      renamed <- c(renamed, sprintf("'%s' -> '%s'", orig[hit], target))
      names(df)[hit] <- target
    }
  }
  if (length(renamed)) {
    cat(sprintf("             renamed in %s: %s\n",
                label, paste(renamed, collapse = ", ")))
  }
  df
}

read_checked <- function(path, label, required_cols) {
  if (!file.exists(path)) {
    xl <- list.files(".", pattern = "\\.xlsx?$", ignore.case = TRUE)
    stop(sprintf(
      paste0("\n  Cannot find the %s file: '%s'\n",
             "  Working directory: %s\n",
             "  Spreadsheets here: %s\n",
             "  Fix: cd to the folder holding the data, or pass paths:\n",
             "    Rscript FigS2_mortality_comparison.R --mortality <file> --cytokine <file>"),
      label, path, getwd(),
      if (length(xl)) paste(xl, collapse = ", ") else "(none)"),
      call. = FALSE)
  }
  df <- as.data.frame(read_excel(path))
  names(df) <- trimws(names(df))
  cat(sprintf("  %-10s %-30s %3d rows | columns: %s\n",
              label, basename(path), nrow(df),
              paste(names(df), collapse = ", ")))
  df <- canonicalize(df, label)

  missing <- setdiff(required_cols, names(df))
  if (length(missing)) {
    stop(sprintf(
      paste0("\n  The %s file is missing required column(s): %s\n",
             "  Columns found: %s\n\n",
             "  Names are matched ignoring case and punctuation, so no\n",
             "  column resembles the one(s) above. Rename the column in the\n",
             "  spreadsheet, or add its spelling to ALIASES."),
      label, paste(missing, collapse = ", "),
      paste(names(df), collapse = ", ")), call. = FALSE)
  }
  df
}

## via character, so factors and text convert by value, not by code
num_col <- function(x) suppressWarnings(as.numeric(as.character(x)))

build_data <- function(mortality_path, cytokine_path) {
  cat("Reading input files\n")
  mort <- read_checked(mortality_path, "mortality",
                       c("Isolate", "Mortality", "SOFA"))
  cyt  <- read_checked(cytokine_path, "cytokine",
                       c("Isolate", "IL_1B", "IL6", "TNF"))

  ## Virulence may sit in either workbook; if absent from both, that
  ## comparison is skipped
  vir_in_mort <- "Virulence" %in% names(mort)
  vir_in_cyt  <- "Virulence" %in% names(cyt)
  if (!vir_in_mort && !vir_in_cyt) {
    warning("no Virulence column in either file; ",
            "the murine virulence comparison will be skipped", call. = FALSE)
    cyt$Virulence <- NA_real_
    vir_in_cyt <- TRUE
  }

  ## stray whitespace otherwise breaks the join silently
  mort$Isolate <- trimws(as.character(mort$Isolate))
  cyt$Isolate  <- trimws(as.character(cyt$Isolate))

  ## "D" in the SOFA column becomes NA (one patient lacks a day-0 score)
  mort$SOFA      <- num_col(mort$SOFA)
  mort$Mortality <- num_col(mort$Mortality)
  if (vir_in_mort) mort$Virulence <- num_col(mort$Virulence)

  ## a patient with no recorded outcome cannot be classified
  n0 <- nrow(mort)
  mort <- mort[!is.na(mort$Mortality), , drop = FALSE]
  if (nrow(mort) < n0) {
    cat(sprintf("  dropped %d row(s) with no recorded Mortality\n",
                n0 - nrow(mort)))
  }

  cyt$IL_1B <- num_col(cyt$IL_1B)
  cyt$IL6   <- num_col(cyt$IL6)
  cyt$TNF   <- num_col(cyt$TNF)
  if (vir_in_cyt) cyt$Virulence <- num_col(cyt$Virulence)

  nonpos <- sum(cyt$IL_1B <= 0 | cyt$IL6 <= 0 | cyt$TNF <= 0, na.rm = TRUE)
  if (nonpos > 0) {
    warning(sprintf("%d replicate(s) have a non-positive concentration; log10 gives NaN/-Inf",
                    nonpos), call. = FALSE)
  }

  cyt$log_IL_1B <- log10(cyt$IL_1B)
  cyt$log_IL6   <- log10(cyt$IL6)
  cyt$log_TNF   <- log10(cyt$TNF)
  cyt$log_TOTAL <- cyt$log_IL_1B + cyt$log_IL6 + cyt$log_TNF

  ## isolate means of the log10 values
  num <- c("log_IL_1B", "log_IL6", "log_TNF", "log_TOTAL")
  nm  <- c("Isolate", "IL_1B", "IL6", "TNF", "TOTAL")
  if (vir_in_cyt) { num <- c(num, "Virulence"); nm <- c(nm, "Virulence_cyt") }
  agg <- aggregate(cyt[num], by = list(Isolate = cyt$Isolate), FUN = mean)
  stopifnot(ncol(agg) == length(nm))
  names(agg) <- nm

  d <- merge(mort, agg, by = "Isolate")

  ## the mortality file (the curated clinical sheet) wins when both carry
  ## virulence; disagreements are named rather than silently reconciled
  if (vir_in_mort && vir_in_cyt) {
    both <- !is.na(d$Virulence) & !is.na(d$Virulence_cyt)
    bad  <- both & d$Virulence != d$Virulence_cyt
    if (any(bad)) {
      warning(sprintf(
        paste0("virulence differs between the two files for %d of %d ",
               "isolates; using the mortality file. Disagreements: %s"),
        sum(bad), sum(both), paste(d$Isolate[bad], collapse = ", ")),
        call. = FALSE)
    }
    d$Virulence_cyt <- NULL
  } else if (vir_in_cyt) {
    names(d)[names(d) == "Virulence_cyt"] <- "Virulence"
  }

  cat(sprintf("  isolates: %d with cytokine data, %d with patient records, %d matched\n",
              nrow(agg), nrow(mort), nrow(d)))
  if (nrow(d) == 0L) {
    stop(paste0("\n  The join produced no rows: isolate names do not match ",
                "between files.\n  Mortality examples: ",
                paste(utils::head(sort(unique(mort$Isolate)), 4), collapse = ", "),
                "\n  Cytokine examples:  ",
                paste(utils::head(sort(unique(cyt$Isolate)), 4), collapse = ", ")),
         call. = FALSE)
  }
  if (nrow(d) != 85L) {
    warning(sprintf("expected 85 matched isolates, got %d", nrow(d)),
            call. = FALSE)
  }
  cat("\n")
  d
}

d <- build_data(f_mort, f_cyto)

## drop any metric with no usable data (e.g. a missing Virulence column)
have <- vapply(VARS, function(v) sum(!is.na(d[[v$col]])) > 0L, logical(1))
if (any(!have)) {
  cat(sprintf("Skipping (no data): %s\n",
              paste(vapply(VARS[!have], `[[`, "", "col"), collapse = ", ")))
  VARS <- VARS[have]
}


## ---- 1. Classify ------------------------------------------------------

d$Group <- factor(ifelse(d$Mortality == 1, "Died", "Survived"),
                  levels = c("Survived", "Died"))

cat(strrep("=", 72), "\n")
cat(sprintf("CLASSIFIED %d isolates by source-patient outcome: %d survived, %d died\n",
            nrow(d), sum(d$Group == "Survived"), sum(d$Group == "Died")))
cat(strrep("=", 72), "\n")

write.csv(d[c("Isolate", "Mortality", "Group", vapply(VARS, `[[`, "", "col"))],
          "classified_isolates.csv", row.names = FALSE)


## ---- 2. Test ----------------------------------------------------------

## rank-biserial r from mid-ranks (fast, for the bootstrap);
## positive = higher values among non-survivors
rank_biserial <- function(surv, died) {
  nd <- length(died); ns <- length(surv)
  rk <- rank(c(died, surv))
  2 * (sum(rk[seq_len(nd)]) - nd * (nd + 1) / 2) / (ns * nd) - 1
}
hodges_lehmann <- function(surv, died) median(as.vector(outer(died, surv, "-")))

boot_ci <- function(surv, died, fn, nboot = NBOOT) {
  out <- numeric(nboot)
  for (i in seq_len(nboot)) {
    out[i] <- fn(sample(surv, length(surv), TRUE), sample(died, length(died), TRUE))
  }
  unname(quantile(out, c(0.025, 0.975), na.rm = TRUE))
}

res <- do.call(rbind, lapply(VARS, function(v) {
  sub  <- d[!is.na(d[[v$col]]), ]
  surv <- as.numeric(sub[[v$col]][sub$Group == "Survived"])
  died <- as.numeric(sub[[v$col]][sub$Group == "Died"])

  p_mw <- suppressWarnings(wilcox.test(surv, died, exact = FALSE)$p.value)
  p_t  <- t.test(surv, died)$p.value                     # Welch, reference only
  ## Shapiro-Wilk on within-group residuals, not the pooled values: a
  ## real group difference would otherwise look like non-normality
  p_sh <- shapiro.test(c(surv - mean(surv), died - mean(died)))$p.value

  r  <- rank_biserial(surv, died);  rci <- boot_ci(surv, died, rank_biserial)
  hl <- hodges_lehmann(surv, died); hci <- boot_ci(surv, died, hodges_lehmann)

  data.frame(
    column = v$col, variable = gsub("\n", " ", v$lab), control = v$ctrl,
    n_surv = length(surv), n_died = length(died),
    median_surv = median(surv),
    iqr_surv = sprintf("%.2f-%.2f", quantile(surv, .25), quantile(surv, .75)),
    median_died = median(died),
    iqr_died = sprintf("%.2f-%.2f", quantile(died, .25), quantile(died, .75)),
    p_mannwhitney = p_mw, p_welch_t = p_t, p_shapiro = p_sh,
    r = r, r_lo = rci[1], r_hi = rci[2],
    hl_shift = hl, hl_lo = hci[1], hl_hi = hci[2],
    stringsAsFactors = FALSE)
}))

## Holm across the laboratory metrics (SOFA held out by default)
res$p_holm <- NA_real_
idx <- if (HOLM_INCLUDES_SOFA) seq_len(nrow(res)) else which(!res$control)
res$p_holm[idx] <- p.adjust(res$p_mannwhitney[idx], method = "holm")

write.csv(res, "mortality_comparison.csv", row.names = FALSE)

fmt_p <- function(p) ifelse(is.na(p), "",
                     ifelse(p < 0.001, "p < 0.001",
                     ifelse(p < 0.1, sprintf("p = %.3f", p), sprintf("p = %.2f", p))))

cat("\nGROUP COMPARISONS (positive r = higher among non-survivors)\n\n")
cat(sprintf("%-28s%16s%16s%12s%12s\n",
            "metric", "survived", "died", "p (MW)", "p (Holm)"))
for (i in seq_len(nrow(res))) {
  x <- res[i, ]
  cat(sprintf("%-28s%9.2f [%s]%9.2f [%s]%12s%12s%s\n",
              x$variable, x$median_surv, x$iqr_surv,
              x$median_died, x$iqr_died,
              fmt_p(x$p_mannwhitney),
              if (is.na(x$p_holm)) "-" else fmt_p(x$p_holm),
              if (x$control) "   (positive control)" else ""))
}
cat("\n")
cat(sprintf("%-28s%26s%24s\n", "metric", "rank-biserial r (95% CI)",
            "median shift (95% CI)"))
for (i in seq_len(nrow(res))) {
  x <- res[i, ]
  cat(sprintf("%-28s%10.3f [%+.3f, %+.3f]%10.3f [%+.3f, %+.3f]\n",
              x$variable, x$r, x$r_lo, x$r_hi, x$hl_shift, x$hl_lo, x$hl_hi))
}
cat("\nWelch t-test p-values (reference only; Shapiro-Wilk p for normality):\n")
for (i in seq_len(nrow(res))) {
  cat(sprintf("  %-28s t: %-12s  shapiro: %.2e%s\n", res$variable[i],
              fmt_p(res$p_welch_t[i]), res$p_shapiro[i],
              if (res$p_shapiro[i] < 0.05) "  (non-normal)" else ""))
}


## ---- Self-check against the Python implementation ---------------------
cat("\n", strrep("=", 72), "\n", sep = "")
cat("SELF-CHECK against the validated Python implementation\n")
cat(strrep("=", 72), "\n")
expected <- data.frame(
  column = c("SOFA", "Virulence", "TOTAL", "IL_1B", "IL6", "TNF"),
  p = c(3.95e-07, 0.1870, 0.3447, 0.2713, 0.6317, 0.3589),
  r = c(0.668, -0.165, -0.124, -0.145, -0.063, -0.121),
  stringsAsFactors = FALSE)
ok <- TRUE
for (i in seq_len(nrow(expected))) {
  got <- res[res$column == expected$column[i], ]
  if (nrow(got) != 1L) {
    cat(sprintf("  %-12s not computed - skipped\n", expected$column[i])); next
  }
  pass <- abs(got$r - expected$r[i]) < 0.005 &&
          abs(log10(got$p_mannwhitney) - log10(expected$p[i])) < 0.02
  ok <- ok && pass
  cat(sprintf("  %-12s p %.3e vs %.3e | r %+.3f vs %+.3f  %s\n",
              expected$column[i], got$p_mannwhitney, expected$p[i],
              got$r, expected$r[i], if (pass) "OK" else "MISMATCH"))
}
cat(if (ok) "\n  All values reproduced.\n" else
    "\n  MISMATCH: check the input files and the virulence source.\n")


## ---- 3. Figure --------------------------------------------------------

make_panel <- function(v, row, tag) {
  sub <- d[!is.na(d[[v$col]]), c("Group", v$col)]
  names(sub)[2] <- "y"

  rng  <- range(sub$y)
  span <- diff(rng)
  if (span == 0) span <- 1
  y_br  <- rng[2] + 0.05 * span      # bracket shoulders
  y_top <- rng[2] + 0.11 * span      # bracket bar
  bracket <- data.frame(x = c(1, 1, 2, 2), y = c(y_br, y_top, y_top, y_br))

  lbl <- sprintf("%s\nn = %d", levels(sub$Group),
                 c(sum(sub$Group == "Survived"), sum(sub$Group == "Died")))
  sig <- !is.na(row$p_mannwhitney) && row$p_mannwhitney < 0.05

  ggplot(sub, aes(Group, y)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.4,
                 colour = "grey25", fill = NA) +
    geom_jitter(aes(colour = Group), width = 0.17, height = 0,
                size = 1.5, alpha = 0.55, show.legend = FALSE) +
    geom_path(data = bracket, aes(x, y), inherit.aes = FALSE,
              colour = "grey40", linewidth = 0.35) +
    annotate("text", x = 1.5, y = y_top + 0.02 * span, vjust = 0,
             label = fmt_p(row$p_mannwhitney), size = 3,
             colour = if (sig) COL_DIED else "grey35",
             fontface = if (sig) "bold" else "plain") +
    scale_colour_manual(values = c(Survived = COL_SURV, Died = COL_DIED)) +
    scale_x_discrete(labels = lbl) +
    expand_limits(y = rng[2] + 0.26 * span) +
    labs(x = NULL, y = v$lab, tag = tag,
         title = if (v$ctrl) "positive control" else NULL) +
    theme_classic(base_size = 9) +
    theme(axis.text  = element_text(size = 8),
          axis.title = element_text(size = 8.5),
          plot.title = element_text(size = 7.5, face = "italic",
                                    colour = "grey45", hjust = 1),
          ## panel letter at the top-left of the whole panel, on both the
          ## combined figure and the individually saved panels
          plot.tag  = element_text(size = 11, face = "bold", hjust = 0, vjust = 1),
          plot.tag.position = "topleft",
          plot.margin = margin(6, 6, 4, 4))
}

## letters run A, B, C ... in VARS order, closing up if a metric was skipped
panels <- lapply(seq_along(VARS),
                 function(i) make_panel(VARS[[i]], res[i, ], LETTERS[i]))

## ragg (or cairo) and cairo_pdf render the Greek glyphs; the stock
## Windows png and pdf devices do not
save_png <- function(path, plot, width, height, dpi = 300) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(path, width = width, height = height, units = "in", res = dpi)
  } else {
    png(path, width = width, height = height, units = "in", res = dpi, type = "cairo")
  }
  on.exit(dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}

dir.create("panels", showWarnings = FALSE)
for (i in seq_along(VARS)) {
  save_png(file.path("panels", paste0(VARS[[i]]$col, ".png")), panels[[i]], 3.2, 3.2)
}

## each panel carries its own letter, so no plot_annotation(tag_levels)
n_col <- if (length(panels) >= 5) 3 else length(panels)
figS2 <- wrap_plots(panels, ncol = n_col)
W <- 3.2 * n_col
H <- 3.2 * ceiling(length(panels) / n_col)

save_png("FigS2_mortality_comparison.png", figS2, W, H)
cairo_pdf("FigS2_mortality_comparison.pdf", width = W, height = H)
print(figS2); invisible(dev.off())

cat("\nwrote classified_isolates.csv, mortality_comparison.csv, panels/,\n",
    "      FigS2_mortality_comparison.png and .pdf\n", sep = "")
