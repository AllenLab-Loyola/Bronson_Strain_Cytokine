# =============================================================================
# LOESS regression of in vitro cytokine pairs (log10-transformed data)
#   IL-1β vs IL-6, IL-1β vs TNF-α, IL-6 vs TNF-α
#
# Workflow:
#   1. Read raw data
#   2. log10-transform every raw replicate value
#   3. (optional) average the log10 replicates per isolate
#   4. LOESS regression + Spearman correlation on the log10 values
#
# Input : cytokine.xlsx  (columns: ID, Isolate, Replicate, IL_1B, IL6, TNF,
#                         ExoS/U, Virulence)
# Output: cytokine_log10_data.csv          (log10-transformed analysis data)
#         cytokine_loess_pairs.pdf / .png  (3-panel figure)
#         cytokine_loess_stats.csv         (Spearman correlation per pair)
# =============================================================================

# ---- Packages ---------------------------------------------------------------
pkgs <- c("readxl", "dplyr", "ggplot2", "patchwork")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

# ---- Settings ---------------------------------------------------------------
infile       <- "cytokine.xlsx"   # path to the raw data file
average_reps <- TRUE    # TRUE: one point per isolate (mean of log10 replicates)
                        # FALSE: plot every replicate as its own point
loess_span   <- 0.75    # LOESS smoothing span (R default = 0.75)
show_ci      <- TRUE    # shaded 95% confidence band around the LOESS curve

cytokines <- c("IL_1B", "IL6", "TNF")

# ---- Read data --------------------------------------------------------------
raw <- read_excel(infile) %>%
  rename(ExoSU = `ExoS/U`) %>%
  mutate(across(all_of(cytokines), as.numeric))

# ---- Step 1: log10-transform the raw data ----------------------------------
# log10 is undefined for values <= 0; flag them rather than silently dropping
n_nonpos <- sum(raw[cytokines] <= 0, na.rm = TRUE)
if (n_nonpos > 0) {
  warning(n_nonpos, " value(s) <= 0 found; set to NA before log10 transform.")
}

dat_log <- raw %>%
  mutate(across(all_of(cytokines), ~ log10(ifelse(.x > 0, .x, NA))))

# ---- Step 2: summarise replicates on the log10 scale ------------------------
# Mean of log10 values = log10 of the geometric mean of the replicates
if (average_reps) {
  dat <- dat_log %>%
    group_by(Isolate, ExoSU, Virulence) %>%
    summarise(across(all_of(cytokines), ~ mean(.x, na.rm = TRUE)),
              .groups = "drop")
} else {
  dat <- dat_log
}

cat("Points per panel:", nrow(dat),
    if (average_reps) "(isolate means of log10 values)" else
      "(individual log10 replicates)", "\n")

write.csv(dat, "cytokine_log10_data.csv", row.names = FALSE)

# ---- Pairs to plot ----------------------------------------------------------
# plotmath expressions render Greek letters on every graphics device/font
labels <- list(
  IL_1B = expression(log[10] ~ "IL-1" * beta ~ "(pg/mL)"),
  IL6   = expression(log[10] ~ "IL-6" ~ "(pg/mL)"),
  TNF   = expression(log[10] ~ "TNF-" * alpha ~ "(pg/mL)")
)

pairs <- list(c("IL_1B", "IL6"),
              c("IL_1B", "TNF"),
              c("IL6",   "TNF"))

# ---- Analysis + plot function ----------------------------------------------
loess_plot <- function(df, x, y) {
  df <- df[complete.cases(df[, c(x, y)]), ]

  # LOESS fit on log10 data (same model geom_smooth draws)
  fit <- loess(reformulate(x, y), data = df, span = loess_span)

  # Spearman correlation on log10 data
  ct <- suppressWarnings(cor.test(df[[x]], df[[y]], method = "spearman"))
  p_txt <- if (ct$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", ct$p.value)
  lab <- bquote(atop("Spearman" ~ rho == .(sprintf("%.2f", ct$estimate)),
                     .(p_txt)))

  g <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(shape = 21, size = 2, fill = "grey70", colour = "grey30",
               alpha = 0.8, stroke = 0.3) +
    geom_smooth(method = "loess", formula = y ~ x, span = loess_span,
                se = show_ci, colour = "#B2182B", fill = "#B2182B",
                alpha = 0.15, linewidth = 1, linetype = "solid") +
    annotation_custom(grid::textGrob(
      as.expression(lab), x = unit(0.04, "npc"), y = unit(0.97, "npc"),
      hjust = 0, vjust = 1, gp = grid::gpar(fontsize = 10))) +
    labs(x = labels[[x]], y = labels[[y]]) +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(colour = "black"),
          aspect.ratio = 1,
          plot.margin = margin(5.5, 18, 5.5, 5.5))

  list(plot = g,
       stats = data.frame(x = x, y = y, n = nrow(df),
                          rho = unname(ct$estimate), p_value = ct$p.value,
                          loess_span = loess_span,
                          loess_residual_se = fit$s))
}

# ---- Build figure -----------------------------------------------------------
res   <- lapply(pairs, function(p) loess_plot(dat, p[1], p[2]))
plots <- lapply(res, `[[`, "plot")
stats <- do.call(rbind, lapply(res, `[[`, "stats"))

fig <- wrap_plots(plots, nrow = 1)

print(fig)
ggsave("cytokine_loess_pairs.pdf", fig, width = 13, height = 4.6,
       device = cairo_pdf)
ggsave("cytokine_loess_pairs.png", fig, width = 13, height = 4.6, dpi = 300)

print(stats)
write.csv(stats, "cytokine_loess_stats.csv", row.names = FALSE)
