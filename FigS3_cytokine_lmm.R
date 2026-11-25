# =============================================================================
# Linear mixed models: log10 cytokine concentration ~ T3SS effector (ExoS vs ExoU)
#
# Data: cytokine.xlsx (Sheet1) - one row per isolate x biological replicate
#   Columns: ID, Isolate, Replicate, IL_1B, IL6, TNF, ExoS/U, Virulence
#
# Model (fit separately for each cytokine):
#   log10(concentration) ~ Effector + (1 | Isolate)
#   - Effector (ExoS / ExoU) is a fixed effect (strain-level trait)
#   - Isolate is a random intercept, accounting for the 3 replicates per
#     strain being non-independent (avoids pseudoreplication)
# =============================================================================

# ---- Packages ---------------------------------------------------------------
pkgs <- c("readxl", "dplyr", "tidyr", "lme4", "lmerTest", "emmeans", "ggplot2")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)

library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)   # Satterthwaite df and p-values for lmer fixed effects
library(emmeans)
library(ggplot2)

# ---- Import -----------------------------------------------------------------
file_path <- "cytokine.xlsx"          # <- change to your file location
raw <- read_excel(file_path, sheet = 1)

dat <- raw %>%
  rename(Effector = `ExoS/U`) %>%
  mutate(
    Isolate   = factor(Isolate),
    Replicate = factor(Replicate),
    Effector  = factor(Effector, levels = c("S", "U"),
                       labels = c("ExoS", "ExoU"))   # ExoS = reference level
  )

# Quick checks
cat("Rows:", nrow(dat), " Isolates:", nlevels(dat$Isolate), "\n")
print(dat %>% distinct(Isolate, Effector) %>% count(Effector, name = "n_isolates"))

# ---- Reshape to long format and log10-transform -----------------------------
cytokine_cols <- c(IL_1B = "IL-1beta", IL6 = "IL-6", TNF = "TNF-alpha")

long <- dat %>%
  pivot_longer(all_of(names(cytokine_cols)),
               names_to = "Cytokine", values_to = "Conc") %>%
  mutate(Cytokine = factor(cytokine_cols[Cytokine], levels = cytokine_cols))

# log10 requires strictly positive values
if (any(long$Conc <= 0, na.rm = TRUE)) {
  stop("Zero/negative concentrations present - decide on a handling rule ",
       "(e.g., substitute LOD/2) before log-transforming.")
}
long <- long %>% mutate(log10Conc = log10(Conc))

# ---- Fit one LMM per cytokine -----------------------------------------------
fit_lmm <- function(d) {
  lmer(log10Conc ~ Effector + (1 | Isolate), data = d, REML = TRUE)
}

models <- long %>%
  split(.$Cytokine) %>%
  lapply(fit_lmm)

# ---- Summaries --------------------------------------------------------------
results <- lapply(names(models), function(cyt) {
  m   <- models[[cyt]]
  cat("\n=====================================================\n")
  cat(" ", cyt, ": log10(conc) ~ Effector + (1 | Isolate)\n")
  cat("=====================================================\n")
  print(summary(m))
  print(anova(m))                                  # Satterthwaite F-test

  # Fixed effect: ExoU - ExoS on log10 scale
  fe  <- summary(m)$coefficients["EffectorExoU", ]
  ci  <- confint(m, parm = "EffectorExoU", method = "profile", quiet = TRUE)

  # Variance components / ICC
  vc  <- as.data.frame(VarCorr(m))
  v_iso <- vc$vcov[vc$grp == "Isolate"]
  v_res <- vc$vcov[vc$grp == "Residual"]

  # Estimated marginal means, back-transformed to concentration (geometric means)
  em  <- emmeans(m, ~ Effector)
  emd <- as.data.frame(em)

  data.frame(
    Cytokine          = cyt,
    GeoMean_ExoS      = 10^emd$emmean[emd$Effector == "ExoS"],
    GeoMean_ExoU      = 10^emd$emmean[emd$Effector == "ExoU"],
    log10_diff_U_vs_S = fe["Estimate"],
    SE                = fe["Std. Error"],
    df                = fe["df"],
    t                 = fe["t value"],
    p_value           = fe["Pr(>|t|)"],
    FoldChange_U_vs_S = 10^fe["Estimate"],
    FC_CI95_low       = 10^ci[1],
    FC_CI95_high      = 10^ci[2],
    Var_Isolate       = v_iso,
    Var_Residual      = v_res,
    ICC               = v_iso / (v_iso + v_res),
    row.names = NULL
  )
})
results <- bind_rows(results) %>%
  mutate(p_adj_Holm = p.adjust(p_value, method = "holm"))  # 3 cytokines tested

cat("\n\n===== Summary: ExoU vs ExoS (reference) =====\n")
print(results %>% mutate(across(where(is.numeric), ~ signif(.x, 4))), row.names = FALSE)
write.csv(results, "LMM_cytokine_effector_results.csv", row.names = FALSE)

# ---- Model diagnostics ------------------------------------------------------
pdf("LMM_diagnostics.pdf", width = 9, height = 6)
for (cyt in names(models)) {
  m <- models[[cyt]]
  par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
  plot(fitted(m), resid(m), xlab = "Fitted", ylab = "Residual", main = "Residuals vs fitted")
  abline(h = 0, lty = 2)
  qqnorm(resid(m), main = "Residual Q-Q"); qqline(resid(m))
  re <- ranef(m)$Isolate[, 1]
  qqnorm(re, main = "Random intercept Q-Q"); qqline(re)
  mtext(cyt, outer = TRUE, font = 2)
}
dev.off()

# ---- Plot: isolate means with model estimates -------------------------------
iso_means <- long %>%
  group_by(Cytokine, Isolate, Effector) %>%
  summarise(log10Conc = mean(log10Conc), .groups = "drop")

emm_df <- bind_rows(lapply(names(models), function(cyt) {
  as.data.frame(emmeans(models[[cyt]], ~ Effector)) %>% mutate(Cytokine = cyt)
})) %>% mutate(Cytokine = factor(Cytokine, levels = cytokine_cols))

p_lab <- results %>%
  mutate(Cytokine = factor(Cytokine, levels = cytokine_cols),
         label = sprintf("p = %.3g", p_value))

# Panel titles with Greek symbols (plotmath expressions, rendered by label_parsed)
greek_labels <- c("IL-1beta"  = 'bold("IL-1")*beta',
                  "IL-6"      = 'bold("IL-6")',
                  "TNF-alpha" = 'bold("TNF-")*alpha')

p <- ggplot(iso_means, aes(Effector, log10Conc)) +
  geom_jitter(colour = "grey25", alpha = 0.7, size = 1.8,
              position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +  # headroom for p-value
  geom_pointrange(data = emm_df,
                  aes(y = emmean, ymin = lower.CL, ymax = upper.CL),
                  size = 0.6, linewidth = 0.9) +
  geom_text(data = p_lab, aes(x = 1.5, y = Inf, label = label),
            vjust = 1.5, inherit.aes = FALSE, size = 3.5) +
  facet_wrap(~ Cytokine, scales = "free_y",
             labeller = as_labeller(greek_labels, default = label_parsed)) +
  labs(x = "T3SS effector", y = expression(log[10]~"concentration (pg/mL)")) +
  theme_bw(base_size = 12) +
  theme(legend.position  = "none",
        strip.background = element_blank(),   # no grey box behind panel titles
        strip.text       = element_text(size = 13),
        panel.grid       = element_blank(),   # no gridlines
        panel.border     = element_blank(),   # remove top/right box lines
        axis.line        = element_line(colour = "black"))  # keep x and y axes

ggsave("LMM_cytokine_effector_plot.pdf", p, width = 9, height = 4)
ggsave("LMM_cytokine_effector_plot.png", p, width = 9, height = 4, dpi = 300)

sessionInfo()
