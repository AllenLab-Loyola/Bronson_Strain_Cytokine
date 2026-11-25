## =====================================================================
## Figure 4: paired clinical analysis -- effect sizes, power, and
## whether day-0 SOFA discriminates mortality better than lab metrics
## Bronson et al., strain-dependent Pa bacteremia manuscript
##
## R translation of clinical_power_analysis.py and compare_effects.py.
##
##   1. Effect sizes: rank-biserial r with bootstrap CIs, AUC,
##      Hodges-Lehmann shifts, odds ratios per SD
##   2. Design-based power: simulated Mann-Whitney power curve and the
##      80% detection limit
##   3. Paired comparisons (DeLong; nested logistic LRT)
##   4. Figure 4: A power curve, B forest plot vs. detection limit
##
## USAGE
##   Rscript Fig4_clinical_power_analysis.R
##   Rscript Fig4_clinical_power_analysis.R --mortality X.xlsx --cytokine Y.xlsx
##
## INPUTS   patient.xlsx (mortality, SOFA, mouse virulence), cytokine.xlsx
## OUTPUTS  effect_sizes_R.csv, power_curve_R.csv, auc_vs_chance_R.csv,
##          delong_comparisons_R.csv, incremental_value_R.csv,
##          Fig4_power_R.png, Fig4_power_R.pdf
##
## DEPENDENCIES  readxl, ggplot2, patchwork, pROC; ragg recommended for
##               Greek glyphs in the PNG (cairo is used otherwise)
##
## RUNTIME
##   About 1-3 minutes. The power simulation dominates (48,000
##   wilcox.test calls). Lower NSIM to 1000 to iterate quickly; raise
##   NBOOT to 50000 to finalise numbers for publication.
##
## REPRODUCIBILITY
##   Bootstrap and power estimates use R's RNG, which differs from the
##   Python implementation's: CIs agree to roughly two decimal places,
##   power to about +/-0.015. Point estimates (r, AUC, HL, OR, p) and
##   the DeLong AUC differences are deterministic; the self-checks
##   verify them against the Python values.
##
## PIPELINE
##   1. log10-transform each replicate measurement
##   2. TOTAL = log10(IL_1B) + log10(IL6) + log10(TNF), per replicate
##   3. mean across the 3 replicates -> one value per isolate
##   4. inner join to patient mortality -> 85 paired cases
##   Rows with no recorded Mortality are dropped before the join.
##   Virulence is read from whichever workbook carries it; if both do,
##   the mortality file wins and disagreements are named in a warning.
##
## TWO DISCREPANCIES AGAINST THE CURRENT DRAFT -- RESOLVE BEFORE SUBMISSION
##   (a) DAY-0 SOFA. The data give medians 4.0 (survivors) vs 9.0
##       (non-survivors), U = 269, p = 4.0e-07, r = 0.67. The draft
##       reports medians 3.0 vs 5.5 and U = 424 (r = 0.48). No subset of
##       these data reproduces the draft values. The data give a
##       STRONGER SOFA effect, but the source of the draft's numbers
##       should be identified.
##   (b) MURINE VIRULENCE SCORE. None of the versions matches the draft
##       (U = 711, p = 0.34):
##         mortality file, "Mouse Virulence":  U = 698, p = 0.19
##         cytokine workbook, "Virulence":     U = 719, p = 0.26
##       The two files disagree for 13 of 85 isolates. This script uses
##       the mortality file. Decide which source is authoritative
##       before reporting that row.
##
## STATISTICAL NOTES
##   * Post-hoc "observed power" is deliberately NOT computed. It is a
##     monotone function of the p-value and adds no information. Use the
##     confidence intervals (what the data exclude) together with the
##     design-based power curve (which depends only on n and alpha).
##   * OR per 1 SD of a standardized continuous predictor corresponds to
##     exp(d), NOT Chinn's exp(d*pi/sqrt(3)); the latter is for a binary
##     exposure contrast. Using it inflates OR per SD from ~1.9 to ~3.2
##     at d = 0.65.
##   * Sign convention: POSITIVE = higher value among non-survivors.
##   * Why a paired comparison: each effect being significant or not
##     against zero does NOT show that two effects differ (Gelman & Stern
##     2006, Am Stat 60:328-331). All predictors are measured on the same
##     patients, so SOFA is compared with each directly using DeLong's
##     test for correlated AUCs (DeLong et al. 1988, Biometrics
##     44:837-845).
##   * pROC: direction = "<" is set explicitly on every roc() call. The
##     default "auto" silently flips any curve with AUC < 0.5, which would
##     report the cytokine AUCs as ~0.56 instead of ~0.44 and invalidate
##     the comparison with SOFA.
## =====================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(pROC)
})


## ---- Settings ---------------------------------------------------------

set.seed(20261118)
NBOOT <- 20000L
NSIM  <- 4000L

## predictors: internal ASCII names (CSV, self-checks) and display labels
VARS <- list(SOFA = "Day-0 SOFA", Virulence = "Murine virulence score",
             TOTAL = "Total cytokine", IL_1B = "IL-1beta",
             IL6 = "IL-6", TNF = "TNFalpha")
REF   <- "SOFA"                                    # reference predictor
GREEK <- c("IL-1beta" = "IL-1β", "TNFalpha" = "TNFα")

## Figure 4 panel B options
##   P_SOURCE: optional p-value column to the right of the forest plot
##     NA            - none (default)
##     "mannwhitney" - does this effect differ from zero? (UNADJUSTED)
##     "delong"      - does this effect differ from SOFA's? (Holm-adjusted)
##   SHOW_EFFECT_COL: print "r [lo, hi]" to the right of each row
P_SOURCE        <- NA
SHOW_EFFECT_COL <- FALSE

## ---- Arguments --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
f_mort <- getarg("--mortality", "patient.xlsx")
f_cyto <- getarg("--cytokine",  "cytokine.xlsx")


## ---- Data -------------------------------------------------------------
## Defensive loader: every step is checked and reported, so a failure
## names the problem rather than surfacing as an opaque error.

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
             "    Rscript Fig4_clinical_power_analysis.R --mortality <file> --cytokine <file>"),
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
  ## analysis is skipped
  vir_in_mort <- "Virulence" %in% names(mort)
  vir_in_cyt  <- "Virulence" %in% names(cyt)
  if (!vir_in_mort && !vir_in_cyt) {
    warning("no Virulence column in either file; ",
            "the murine virulence analysis will be skipped", call. = FALSE)
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

  ## patients with no recorded outcome cannot contribute to any analysis
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


## ---- Effect-size and power functions ----------------------------------

## Rank-biserial correlation from mid-ranks: ~100x faster than
## wilcox.test, which matters inside the bootstrap. Identical to
## 2 * W / (n1 * n2) - 1 with W from wilcox.test(died, surv), ties
## included. Positive = higher values among non-survivors.
rank_biserial <- function(surv, died) {
  nd <- length(died); ns <- length(surv)
  rk <- rank(c(died, surv))                      # average ranks for ties
  u  <- sum(rk[seq_len(nd)]) - nd * (nd + 1) / 2 # Mann-Whitney U for died
  2 * u / (ns * nd) - 1
}

hodges_lehmann <- function(surv, died) {
  median(as.vector(outer(died, surv, "-")))
}

boot_ci <- function(surv, died, fn, nboot = NBOOT) {
  out <- numeric(nboot)
  ns <- length(surv); nd <- length(died)
  for (i in seq_len(nboot)) {
    out[i] <- fn(sample(surv, ns, replace = TRUE),
                 sample(died, nd, replace = TRUE))
  }
  unname(quantile(out, c(0.025, 0.975), na.rm = TRUE))
}

effect_table <- function(d, vars) {
  res <- do.call(rbind, lapply(names(vars), function(v) {
    sub  <- d[!is.na(d[[v]]), ]
    surv <- as.numeric(sub[[v]][sub$Mortality == 0])
    died <- as.numeric(sub[[v]][sub$Mortality == 1])

    p_mw <- suppressWarnings(wilcox.test(surv, died, exact = FALSE)$p.value)
    r    <- rank_biserial(surv, died)

    ## confirm the fast computation agrees with wilcox.test (once per
    ## variable, outside the bootstrap)
    w_chk <- as.numeric(suppressWarnings(
      wilcox.test(died, surv, exact = FALSE)$statistic))
    r_chk <- 2 * w_chk / (length(surv) * length(died)) - 1
    if (abs(r - r_chk) > 1e-9) {
      stop(sprintf("rank_biserial disagrees with wilcox.test for %s: %.6f vs %.6f",
                   vars[[v]], r, r_chk))
    }

    r_ci  <- boot_ci(surv, died, rank_biserial)
    hl    <- hodges_lehmann(surv, died)
    hl_ci <- boot_ci(surv, died, hodges_lehmann)

    x   <- as.vector(scale(as.numeric(sub[[v]])))
    co  <- summary(glm(sub$Mortality ~ x, family = binomial()))$coefficients
    b   <- co[2, 1]; se <- co[2, 2]

    data.frame(
      variable = vars[[v]],
      n = length(surv) + length(died),
      n_surv = length(surv), n_died = length(died),
      median_surv = median(surv), median_died = median(died),
      p_mannwhitney = p_mw,
      r = r, r_lo = r_ci[1], r_hi = r_ci[2],
      auc = (r + 1) / 2, auc_lo = (r_ci[1] + 1) / 2,
      auc_hi = (r_ci[2] + 1) / 2,
      hl_shift = hl, hl_lo = hl_ci[1], hl_hi = hl_ci[2],
      or_per_sd = exp(b), or_lo = exp(b - 1.96 * se),
      or_hi = exp(b + 1.96 * se), p_logistic = co[2, 4],
      stringsAsFactors = FALSE)
  }))
  rownames(res) <- NULL
  res
}

power_mannwhitney <- function(d_true, n_surv, n_died,
                              nsim = NSIM, alpha = 0.05) {
  hits <- 0L
  for (i in seq_len(nsim)) {
    s <- rnorm(n_surv)
    k <- rnorm(n_died, mean = d_true)
    p <- suppressWarnings(wilcox.test(s, k, exact = FALSE)$p.value)
    if (!is.na(p) && p < alpha) hits <- hits + 1L
  }
  hits / nsim
}

power_logistic <- function(or_sd, n, event_rate,
                           nsim = NSIM, alpha = 0.05) {
  b1 <- log(or_sd)
  z  <- rnorm(200000)
  f  <- function(b0) mean(plogis(b0 + b1 * z)) - event_rate
  b0 <- uniroot(f, c(-10, 10))$root
  hits <- 0L
  for (i in seq_len(nsim)) {
    x <- rnorm(n)
    y <- rbinom(n, 1, plogis(b0 + b1 * x))
    if (sum(y) > 0 && sum(y) < n) {
      fit <- try(suppressWarnings(glm(y ~ x, family = binomial())),
                 silent = TRUE)
      if (!inherits(fit, "try-error")) {
        p <- summary(fit)$coefficients[2, 4]
        if (!is.na(p) && p < alpha) hits <- hits + 1L
      }
    }
  }
  hits / nsim
}

d_to_auc <- function(d) pnorm(d / sqrt(2))
d_to_r   <- function(d) 2 * d_to_auc(d) - 1

## n per group for a two-sample t-test at given power (Cohen's d)
n_for_d <- function(d, ratio, power = 0.80, alpha = 0.05) {
  f <- function(n1) {
    n2 <- n1 * ratio
    ncp <- d / sqrt(1 / n1 + 1 / n2)
    crit <- qt(1 - alpha / 2, n1 + n2 - 2)
    (pt(-crit, n1 + n2 - 2, ncp) + pt(crit, n1 + n2 - 2, ncp,
                                      lower.tail = FALSE)) - power
  }
  uniroot(f, c(3, 100000))$root
}

rule <- function(w = 78) cat(strrep("=", w), "\n", sep = "")


## ---- 1. Data and observed effects -------------------------------------
## The RNG-consuming steps (bootstrap, then power curve, then logistic
## cross-check) run in this order; changing the order changes the CIs
## and power estimates.

d <- build_data(f_mort, f_cyto)

## drop any variable with no usable data (e.g. a missing Virulence column)
usable <- vapply(names(VARS), function(v) sum(!is.na(d[[v]])) > 0L, logical(1))
if (any(!usable)) {
  cat(sprintf("Skipping (no data): %s\n\n", paste(unlist(VARS[!usable]), collapse = ", ")))
  VARS <- VARS[usable]
}

n_died <- sum(d$Mortality == 1)
n_surv <- sum(d$Mortality == 0)

rule()
cat(sprintf("PAIRED COHORT: n = %d (%d died, %d survived); SOFA for %d\n",
            nrow(d), n_died, n_surv, sum(!is.na(d$SOFA))))
rule(); cat("\n")

res <- effect_table(d, VARS)

cat("OBSERVED EFFECTS (positive = higher in non-survivors)\n\n")
for (i in seq_len(nrow(res))) {
  x <- res[i, ]
  cat(sprintf("%-24s %5.2f/%-6.2f p=%8.1e  r=%+.3f [%+.3f,%+.3f]  AUC=%.3f\n",
              x$variable, x$median_surv, x$median_died, x$p_mannwhitney,
              x$r, x$r_lo, x$r_hi, x$auc))
}
cat("\n")
for (i in seq_len(nrow(res))) {
  x <- res[i, ]
  cat(sprintf("%-24s HL=%+7.3f [%+.3f,%+.3f]   OR/SD=%5.2f [%.2f,%.2f]\n",
              x$variable, x$hl_shift, x$hl_lo, x$hl_hi,
              x$or_per_sd, x$or_lo, x$or_hi))
}


## ---- 2. Design-based power --------------------------------------------

dgrid <- c(0.1, 0.2, 0.3, 0.35, 0.4, 0.5, 0.6, 0.65, 0.7, 0.8, 0.9, 1.0)
pc <- data.frame(d = dgrid, auc = d_to_auc(dgrid),
                 r = d_to_r(dgrid), or_per_sd = exp(dgrid),
                 power = vapply(dgrid, power_mannwhitney, numeric(1),
                                n_surv = n_surv, n_died = n_died))

## 80% detection limit, interpolated off the simulated curve
mde_d <- approx(pc$power, pc$d, xout = 0.80)$y
mde_r <- d_to_r(mde_d)

cat("\n"); rule()
cat("DESIGN-BASED POWER (depends on n and alpha only, not on the data)\n")
rule()
cat(sprintf("%6s %7s %7s %7s %7s\n", "d", "|r|", "AUC", "OR/SD", "power"))
for (i in seq_len(nrow(pc))) {
  cat(sprintf("%6.2f %7.3f %7.3f %7.2f %7.3f\n",
              pc$d[i], pc$r[i], pc$auc[i], pc$or_per_sd[i], pc$power[i]))
}
cat(sprintf("\n  80%% power requires d = %.2f (|r| = %.2f, AUC = %.2f, OR/SD = %.2f)\n",
            mde_d, mde_r, d_to_auc(mde_d), exp(mde_d)))

ratio <- n_surv / n_died
cat("\n  Paired cases required for 80% power at smaller true effects:\n")
for (dv in c(0.20, 0.30, 0.35)) {
  n1 <- n_for_d(dv, ratio)
  cat(sprintf("    d = %.2f (|r| = %.2f): %5.0f cases (%.0f deaths)\n",
              dv, d_to_r(dv), n1 * (1 + ratio), n1))
}

cat(sprintf("\n  Logistic cross-check, power at OR/SD = 1.94: %.3f\n",
            power_logistic(1.94, nrow(d), n_died / nrow(d))))

sofa     <- res[res$variable == VARS[[REF]], ]
lab_rows <- res[res$variable != VARS[[REF]], ]
cat("\n"); rule()
cat("SUMMARY FOR THE MANUSCRIPT\n")
rule()
cat(sprintf("  Detection limit (80%% power):   |r| = %.2f, OR per SD = %.2f\n",
            mde_r, exp(mde_d)))
cat(sprintf("  Day-0 SOFA:                    |r| = %.2f -> %.1fx the limit (ABOVE)\n",
            abs(sofa$r), abs(sofa$r) / mde_r))
cat(sprintf("  Laboratory metrics:            |r| = %.2f to %.2f -> %.2fx or less (BELOW)\n",
            min(abs(lab_rows$r)), max(abs(lab_rows$r)), max(abs(lab_rows$r)) / mde_r))
cat(sprintf("  Widest laboratory OR/SD CI:    [%.2f, %.2f]\n",
            min(lab_rows$or_lo), max(lab_rows$or_hi)))


## ---- 3. Paired comparisons with SOFA ----------------------------------

## paired comparisons need complete data on every predictor
cc <- d[stats::complete.cases(d[names(VARS)]), ]
y  <- cc$Mortality

cat("\n"); rule(76)
cat(sprintf("Complete-case cohort for paired tests: n = %d (%d died, %d survived)\n",
            nrow(cc), sum(y == 1), sum(y == 0)))
cat(sprintf("(dropped %d of %d for missing data)\n", nrow(d) - nrow(cc), nrow(d)))
rule(76)

## one ROC per predictor, orientation fixed (see header)
rocs <- lapply(names(VARS), function(v)
  roc(response = y, predictor = cc[[v]],
      levels = c(0, 1), direction = "<", quiet = TRUE))
names(rocs) <- names(VARS)
others <- setdiff(names(VARS), REF)

## Q1: does each predictor beat chance? AUC vs 0.5, DeLong variance
q1 <- do.call(rbind, lapply(names(VARS), function(v) {
  ci  <- as.numeric(ci.auc(rocs[[v]], method = "delong"))   # lo, auc, hi
  auc <- as.numeric(rocs[[v]]$auc)
  z   <- (auc - 0.5) / ((ci[3] - ci[1]) / (2 * 1.96))
  data.frame(predictor = VARS[[v]], auc = auc, lo = ci[1], hi = ci[3],
             z = z, p = 2 * pnorm(-abs(z)))
}))

cat("\nQ1. Does each predictor discriminate better than chance?\n")
cat("    AUC vs 0.50, DeLong variance\n\n")
cat(sprintf("%-24s%8s%20s%11s\n", "predictor", "AUC", "95% CI", "p"))
for (i in seq_len(nrow(q1))) {
  cat(sprintf("%-24s%8.3f%20s%11.2e\n", q1$predictor[i], q1$auc[i],
              sprintf("[%.3f, %.3f]", q1$lo[i], q1$hi[i]), q1$p[i]))
}

## Q2 (key test): does SOFA beat each other predictor? paired DeLong
q2 <- do.call(rbind, lapply(others, function(v) {
  tt <- roc.test(rocs[[REF]], rocs[[v]], method = "delong", paired = TRUE)
  data.frame(comparison = VARS[[v]],
             auc_sofa  = as.numeric(rocs[[REF]]$auc),
             auc_other = as.numeric(rocs[[v]]$auc),
             diff = as.numeric(rocs[[REF]]$auc) - as.numeric(rocs[[v]]$auc),
             z = unname(tt$statistic), p = tt$p.value)
}))
q2$p_holm <- p.adjust(q2$p, method = "holm")

cat("\n"); rule(76)
cat("Q2. Does day-0 SOFA discriminate BETTER than each other predictor?\n")
cat("    Paired DeLong test for correlated AUCs, Holm-adjusted\n\n")
cat(sprintf("%-24s%11s%12s%8s%7s%11s%11s\n", "SOFA vs.", "AUC(SOFA)",
            "AUC(other)", "diff", "z", "p", "p(Holm)"))
for (i in seq_len(nrow(q2))) {
  x <- q2[i, ]
  cat(sprintf("%-24s%11.3f%12.3f%8.3f%7.2f%11.2e%11.2e\n",
              x$comparison, x$auc_sofa, x$auc_other, x$diff,
              x$z, x$p, x$p_holm))
}

## Q3: does a predictor add to SOFA? nested logistic models, LRT
z0   <- as.vector(scale(cc[[REF]]))
base <- glm(y ~ z0, family = binomial())
q3 <- do.call(rbind, lapply(others, function(v) {
  x    <- as.vector(scale(cc[[v]]))
  full <- glm(y ~ z0 + x, family = binomial())
  lr   <- as.numeric(2 * (logLik(full) - logLik(base)))
  data.frame(added = VARS[[v]],
             or_per_sd = exp(coef(full)[["x"]]),
             lr_chisq = lr, p = pchisq(lr, df = 1, lower.tail = FALSE))
}))
q3$p_holm <- p.adjust(q3$p, method = "holm")

cat("\n"); rule(76)
cat("Q3. Does any predictor add information to day-0 SOFA?\n")
cat("    Nested logistic models, likelihood-ratio test, Holm-adjusted\n\n")
cat(sprintf("%-24s%11s%10s%9s%10s\n", "added to SOFA", "OR per SD",
            "LR chi2", "p", "p(Holm)"))
for (i in seq_len(nrow(q3))) {
  x <- q3[i, ]
  cat(sprintf("%-24s%11.2f%10.2f%9.3f%10.3f\n",
              x$added, x$or_per_sd, x$lr_chisq, x$p, x$p_holm))
}


## ---- Self-checks against the Python implementation --------------------
## Deterministic values only; bootstrap CIs and power are stochastic.
## Virulence is the "Mouse Virulence" column of the mortality file. The
## cytokine-workbook copy differs for 13 isolates and gives r = -0.141,
## OR = 0.75 and a SOFA AUC difference of 0.395 instead; a mismatch on
## that row alone most likely means the other source was used.
cat("\n"); rule()
cat("SELF-CHECK against the Python implementation (deterministic values)\n")
rule()

expected <- data.frame(
  variable = c("Day-0 SOFA", "Murine virulence score", "Total cytokine",
               "IL-1beta", "IL-6", "TNFalpha"),
  r    = c(0.668, -0.165, -0.124, -0.145, -0.063, -0.121),
  or   = c(5.50, 0.68, 0.88, 0.82, 0.90, 0.99),
  diff = c(NA, 0.407, 0.395, 0.406, 0.364, 0.393),   # SOFA AUC minus this
  stringsAsFactors = FALSE)

ok <- TRUE
for (i in seq_len(nrow(expected))) {
  e   <- expected[i, ]
  got <- res[res$variable == e$variable, ]
  if (nrow(got) != 1L) {
    cat(sprintf("  %-24s not computed - skipped\n", e$variable))
    next
  }
  pass <- abs(got$r - e$r) < 0.005 && abs(got$or_per_sd - e$or) < 0.02
  line <- sprintf("  %-24s r %+.3f vs %+.3f | OR %5.2f vs %5.2f",
                  e$variable, got$r, e$r, got$or_per_sd, e$or)
  if (!is.na(e$diff)) {
    g2   <- q2$diff[q2$comparison == e$variable]
    pass <- pass && length(g2) == 1L && abs(g2 - e$diff) < 0.005
    line <- paste0(line, sprintf(" | AUC diff %.3f vs %.3f", g2, e$diff))
  }
  ok <- ok && pass
  cat(line, if (pass) " OK\n" else " MISMATCH\n")
}
cat(if (ok) "\n  All deterministic values reproduced.\n" else
    paste0("\n  MISMATCH: check the input files, the pipeline, and that ",
           "direction = \"<\" is set on every roc() call.\n"))


## ---- Write results ----------------------------------------------------

write.csv(res, "effect_sizes_R.csv",       row.names = FALSE)
write.csv(pc,  "power_curve_R.csv",        row.names = FALSE)
write.csv(q1,  "auc_vs_chance_R.csv",      row.names = FALSE)
write.csv(q2,  "delong_comparisons_R.csv", row.names = FALSE)
write.csv(q3,  "incremental_value_R.csv",  row.names = FALSE)


## ---- 4. Figure 4 ------------------------------------------------------

RED  <- "#c0392b"   # SOFA, detection limit
BLUE <- "#2c6fbb"   # laboratory metrics
GREY <- "#9a9a9a"   # shading and guide lines
SHADE_ALPHA <- 0.14

## display labels (Greek); rows drawn top to bottom in table order
fig <- res
fig$label <- ifelse(fig$variable %in% names(GREEK),
                    GREEK[fig$variable], fig$variable)
fig$label <- factor(fig$label, levels = rev(fig$label))
fig$grp   <- ifelse(fig$variable == VARS[[REF]], "sofa", "lab")
fig$txt   <- sprintf("%+.2f [%+.2f, %+.2f]", fig$r, fig$r_lo, fig$r_hi)

## optional p-value column (see P_SOURCE in Settings)
fmt_p <- function(p) ifelse(is.na(p), "(reference)",
                     ifelse(p < 0.001, "p < 0.001",
                     ifelse(p < 0.1, sprintf("p = %.3f", p),
                            sprintf("p = %.2f", p))))
SHOW_P <- !is.na(P_SOURCE)
if (!SHOW_P) {
  fig$ptxt <- ""; P_HEAD <- ""
} else if (P_SOURCE == "delong") {
  ## NA on the SOFA row -> "(reference)"
  fig$ptxt <- fmt_p(q2$p_holm[match(fig$variable, q2$comparison)])
  P_HEAD   <- "p vs. SOFA"
} else {
  fig$ptxt <- fmt_p(fig$p_mannwhitney)
  P_HEAD   <- "p vs. 0"
}

## axis.text.y takes a vector in level order (bottom row first);
## only SOFA is highlighted
is_sofa   <- levels(fig$label) == VARS[[REF]]
tick_col  <- ifelse(is_sofa, RED, "black")
tick_face <- ifelse(is_sofa, "bold", "plain")

## Tick labels. Breaks are literal vectors because seq() leaves a
## floating-point residue (~-1e-16) instead of 0; zapsmall() flattens
## any that remains before formatting.
fmt1 <- function(x) ifelse(is.na(x), "", sprintf("%.1f", zapsmall(x, digits = 10)))
BREAKS_A <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
BREAKS_B <- c(-0.5, 0, 0.5, 1.0)

## optional columns sit right of panel B (which spans -0.75 to 1.0);
## RIGHT_PAD (points) reserves room for the longest label drawn
EFFECT_X <- 1.06
P_X      <- if (SHOW_EFFECT_COL) 1.95 else 1.06
pad_chars <- 0
if (SHOW_EFFECT_COL) pad_chars <- pad_chars + max(nchar(fig$txt)) + 2
if (SHOW_P)          pad_chars <- pad_chars + max(nchar(fig$ptxt)) + 4
RIGHT_PAD <- if (pad_chars > 0) pad_chars * 5.6 else 8

## plot.title.position = "plot" aligns titles with the left edge of the
## whole plot, including the axis labels
thm <- theme_classic(base_size = 9) +
  theme(plot.title          = element_text(face = "bold", size = 9.5,
                                           hjust = 0, margin = margin(b = 6)),
        plot.title.position = "plot",
        plot.tag            = element_text(face = "bold", size = 16),
        plot.tag.position   = "topleft",
        axis.text           = element_text(size = 8),
        axis.title          = element_text(size = 9),
        axis.title.x        = element_text(lineheight = 1.15),
        plot.margin         = margin(5, 6, 5, 6))

## panel A: power curve
pA <- ggplot(pc, aes(r, power * 100)) +
  annotate("rect", xmin = 0, xmax = mde_r, ymin = 0, ymax = 103,
           fill = GREY, alpha = SHADE_ALPHA) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = GREY, linewidth = 0.4) +
  geom_vline(xintercept = mde_r, linetype = "dotted", colour = RED, linewidth = 0.6) +
  geom_line(colour = "#222222", linewidth = 0.9) +
  annotate("text", x = 0.02, y = 83, hjust = 0, size = 2.8,
           colour = "#777777", label = "80% power") +
  annotate("text", x = mde_r / 2, y = 97, hjust = 0.5, size = 2.8,
           colour = "#777777", fontface = "italic", label = "underpowered") +
  annotate("text", x = mde_r + 0.015, y = 4, hjust = 0, vjust = 0,
           size = 2.8, colour = RED, lineheight = 0.95,
           label = sprintf("80%% power at\n|r| = %.2f", mde_r)) +
  scale_x_continuous(breaks = BREAKS_A, labels = fmt1) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  ## train the scales over the full displayed range so end breaks show
  expand_limits(x = c(0, 0.55), y = c(0, 103)) +
  coord_cartesian(xlim = c(0, 0.55), ylim = c(0, 103), expand = FALSE) +
  labs(x = "True effect size (|rank-biserial r|)", y = "Power (%)",
       title = "Power of this cohort", tag = "A") +
  thm

## panel B: forest plot; clip = "off" lets the optional columns sit
## outside the panel
pB <- ggplot(fig, aes(r, label, colour = grp)) +
  scale_y_discrete() +      # declare y discrete before the numeric shading rectangle
  annotate("rect", xmin = -mde_r, xmax = mde_r,
           ymin = 0.4, ymax = nrow(fig) + 0.6, fill = GREY, alpha = SHADE_ALPHA) +
  geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.35) +
  geom_vline(xintercept = c(-mde_r, mde_r), linetype = "dotted",
             colour = RED, linewidth = 0.6) +
  geom_errorbar(aes(xmin = r_lo, xmax = r_hi), width = 0.16,
                linewidth = 0.55, orientation = "y") +
  geom_point(size = 2.1) +
  annotate("text", x = 0, y = nrow(fig) + 0.5, hjust = 0.5, size = 2.6,
           colour = "#777777", fontface = "italic",
           label = "shaded = undetectable at 80% power") +
  scale_colour_manual(values = c(sofa = RED, lab = BLUE), guide = "none") +
  scale_x_continuous(breaks = BREAKS_B, labels = fmt1) +
  expand_limits(x = c(-0.75, 1.0)) +
  coord_cartesian(xlim = c(-0.75, 1.0), ylim = c(0.2, nrow(fig) + 0.6),
                  expand = FALSE, clip = "off") +
  labs(x = paste0("Observed effect size (rank-biserial r, 95% CI)\n",
                  "← lower in non-survivors      ",
                  "higher in non-survivors →"),
       y = NULL, title = "Observed effects vs. detection limit", tag = "B") +
  thm +
  theme(axis.line.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(colour = tick_col, face = tick_face, size = 8.5),
        plot.margin  = margin(5, RIGHT_PAD, 5, 6))

if (SHOW_EFFECT_COL) {
  pB <- pB + geom_text(aes(x = EFFECT_X, label = txt), hjust = 0, size = 2.7,
                       show.legend = FALSE)
}
if (SHOW_P) {
  pB <- pB +
    geom_text(aes(x = P_X, label = ptxt), hjust = 0, size = 2.7, show.legend = FALSE) +
    annotate("text", x = P_X, y = nrow(fig) + 0.5, hjust = 0, size = 2.5,
             colour = "#777777", fontface = "italic", label = P_HEAD)
}

## No overall title: the cohort size belongs in the figure legend, e.g.
##   "Resolution of the paired clinical cohort (n = 85; 31 deaths,
##    54 survivors)."
fig4 <- (pA | pB) + plot_layout(widths = c(1, 1.25))

## Devices opened explicitly. ragg (or cairo) and cairo_pdf render the
## Greek glyphs; the stock Windows png and pdf devices do not.
W <- 10.2; H <- 4.0
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_png("Fig4_power_R.png", width = W, height = H, units = "in", res = 300)
} else {
  png("Fig4_power_R.png", width = W, height = H, units = "in", res = 300, type = "cairo")
}
print(fig4); invisible(dev.off())

cairo_pdf("Fig4_power_R.pdf", width = W, height = H)
print(fig4); invisible(dev.off())

cat("\nwrote effect_sizes_R.csv, power_curve_R.csv, auc_vs_chance_R.csv,\n",
    "      delong_comparisons_R.csv, incremental_value_R.csv,\n",
    "      Fig4_power_R.png, Fig4_power_R.pdf\n", sep = "")
