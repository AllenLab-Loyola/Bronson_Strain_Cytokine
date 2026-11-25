## =====================================================================
## Figure 3: cytokine stimulation, virulence and mouse survival
##
##   A  Isolate cytokine level vs. virulence score (linear mixed models)
##   B  Hazard ratio per log10 cytokine, within each dose quartile (Cox)
##   C  Kaplan-Meier survival, low- vs. high-stimulating strains,
##      6.6-7.0 log10 CFU dose quartile
##
## Also writes the supplementary Kaplan-Meier figure (all dose quartiles)
## and the mortality table.
##
## Inputs : cytokine.xlsx, mouse.xlsx
## Outputs: CombinedFig.pdf             Figure 3
##          aim2_km_dose_quartiles.pdf  supplementary KM figure
##          aim2_km_table.csv           mortality by stimulation group
##
## Also needs (called with ::, not attached): car, pbkrtest
## (Kenward-Roger F tests), scales
## =====================================================================

library(tidyverse)
library(readxl)
library(lmerTest)    # lmer() with Kenward-Roger tests; loads lme4
library(emmeans)
library(survival)
library(patchwork)

set.seed(1)          # the jittered isolate points in panel A


## ---- Settings ---------------------------------------------------------

#data-quality outlier, excluded from all analyses
EXCLUDE_ISOLATES = c("PABL99")

#log10 CFU dose bins, roughly quartiles
DOSE_BREAKS = c(5.1, 6.6, 7, 7.4, 9.2)

#dose quartile shown in panel C. Chosen on design grounds, not on
#significance: its observed mortality (49%) is closest to 50%, where a
#two-group survival comparison has the most power, and it brackets the
#nominal ~7.0 log10 CFU challenge dose described in the methods.
FOCUS_Q = "6.6-7.0"

#one row per cytokine: raw column, isolate-mean column, display labels
CYTO = tibble(
  raw  = c("IL_1B",      "IL6",      "TNF"),
  avg  = c("avg_il1b",   "avg_il6",  "avg_tnf"),
  name = c("IL-1β", "IL-6",     "TNF-α"),     # plain text
  expr = c("'IL-1'*beta", "'IL-6'",  "'TNF-'*alpha")    # plotmath
)

#Kaplan-Meier group colours
col_low  = "#1A4683"   # weak stimulators   (matches Fig 1C "Uniform Low")
col_high = "#8B2211"   # strong stimulators (matches Fig 1C "Uniform High")
grp_lab  = c(Low = "Low-Stimulating", High = "High-Stimulating")


## ---- Figure style -----------------------------------------------------
#Change these numbers once and all nine panels follow.
FIG_BASE  = 14                          # axis and panel titles (pt)
FIG_TICK  = 12                          # tick labels, legend (pt)
FIG_ANNOT = 11                          # in-panel text (pt)
ANNOT_SZ  = FIG_ANNOT / ggplot2::.pt    # geom_text size is in mm

theme_fig = function() {
  theme_classic(base_size = FIG_BASE) +
    theme(
      text              = element_text(colour = "black"),
      axis.title        = element_text(size = FIG_BASE, colour = "black"),
      axis.title.x      = element_text(margin = margin(t = 8)),
      axis.title.y      = element_text(margin = margin(r = 8)),
      axis.text         = element_text(size = FIG_TICK, colour = "black"),
      axis.line         = element_line(colour = "black", linewidth = 0.5),
      axis.ticks        = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(4, "pt"),
      plot.title        = element_text(size = FIG_BASE, colour = "black"),
      legend.text       = element_text(size = FIG_TICK, colour = "black"),
      legend.position   = "bottom",
      plot.tag          = element_text(size = 32, face = "bold", colour = "black"),
      plot.margin       = margin(8, 10, 8, 8)
    )
}

#lighter style for the supplementary all-quartile KM figure
theme_journal = function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title       = element_text(colour = "grey10"),
      axis.title.y     = element_text(margin = margin(r = 14)),
      axis.text        = element_text(colour = "grey30"),
      axis.line        = element_line(colour = "grey55", linewidth = 0.3),
      axis.ticks       = element_line(colour = "grey55", linewidth = 0.3),
      panel.grid.major = element_blank(),
      plot.margin      = margin(8, 10, 6, 8)
    )
}

#p-value formatters
fmt_p = function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
fmt_p_expr = function(p)                       # plotmath, for KM panels
  ifelse(p < 0.001, "italic(p)[adj]*' < 0.001'",
         sprintf("italic(p)[adj]*' = %.3f'", p))
fmt_p_tab = function(p) {                      # journal style, for the table
  if (is.na(p)) return(NA_character_)
  if (p >= 0.001) return(sprintf("%.3f", p))
  e = floor(log10(p))
  sprintf("%.2f x 10^%d", p / 10^e, e)
}

#cut() writes bins as [5.1,6.6); relabel as 5.1-6.6
pretty_bins = function(x) {
  s = as.character(x); out = s; iv = grepl(",", s)
  b = gsub("[][()]", "", s[iv])
  out[iv] = vapply(strsplit(b, ","),
                   function(z) sprintf("%.1f-%.1f", as.numeric(z[1]), as.numeric(z[2])),
                   character(1))
  out
}


## ---- Data -------------------------------------------------------------

#replicate-level cytokines, log10 transformed
pa_cyto = read_xlsx("cytokine.xlsx") %>%
  filter(!Isolate %in% EXCLUDE_ISOLATES) %>%
  mutate(across(all_of(CYTO$raw), log10))

#isolate-level means of the log10 values
pa_cyto_iso = pa_cyto %>%
  group_by(ID, Isolate, `ExoS/U`, Virulence) %>%
  summarise(avg_il1b = mean(IL_1B),
            avg_il6  = mean(IL6),
            avg_tnf  = mean(TNF), .groups = "drop")

#mouse experiments: alive counts at each time point -> one row per mouse
#with its death (or censoring) time
expand_mice = function(df, key) {
  final_time = max(df$Time)
  n_mice     = max(df$AliveCount, na.rm = TRUE)
  death_time = map_int(seq_len(n_mice), function(m) {
    idx = which(df$AliveCount < m)
    if (length(idx) == 0) final_time else df$Time[min(idx)]
  })
  tibble(MouseID = seq_len(n_mice),
         Time    = death_time,
         Dead    = as.numeric(death_time < final_time))
}

mouse_long = read_xlsx("mouse.xlsx") %>%
  filter(!is.na(time88), Strain != "BL058") %>%
  pivot_longer(starts_with("time"), names_to = "Time", values_to = "AliveCount") %>%
  mutate(Time = as.numeric(str_remove(Time, "time"))) %>%
  arrange(ID, Strain, Time) %>%
  group_by(ID, Strain, Date, `Log CFU`) %>%
  group_modify(expand_mice) %>%
  ungroup()

cat(sprintf("mouse experiments without cytokine data (dropped): %d\n",
            n_distinct(anti_join(mouse_long, pa_cyto_iso, by = "ID")$ID)))

mouse_mod = inner_join(mouse_long, pa_cyto_iso, by = "ID") %>%
  mutate(
    #dose quartiles, used in panels B and C
    cat_dose2 = cut(`Log CFU`, breaks = DOSE_BREAKS, include.lowest = TRUE,
                    right = FALSE, dig.lab = 4),
    #the three nominal dose targets (Table 1)
    cat_dose  = case_when(`Log CFU` < 6.25 ~ NA_real_,
                          `Log CFU` < 6.75 ~ 6.5,
                          `Log CFU` < 7.25 ~ 7,
                          `Log CFU` < 7.75 ~ 7.75,
                          TRUE             ~ NA_real_)
  )


## ---- Cohort summary ---------------------------------------------------

model_dat = mouse_mod %>% filter(!is.na(cat_dose2))

model_dat %>%
  summarise(n_isolates = n_distinct(ID), n_strains = n_distinct(Strain),
            n_mice = n(), n_deaths = sum(Dead))

model_dat %>%
  group_by(cat_dose2) %>%
  summarise(n_isolates = n_distinct(ID), n_mice = n(), n_deaths = sum(Dead),
            mortality = mean(Dead))

#mice per isolate (cluster balance)
model_dat %>%
  count(ID, name = "mice") %>%
  summarise(min = min(mice), median = median(mice), max = max(mice))


## ---- Panel A: cytokine vs. virulence score ----------------------------

trend_panel = function(i) {
  raw = CYTO$raw[i]; avg = CYTO$avg[i]

  #linear trend: virulence score as a continuous predictor
  lin   = lmer(as.formula(sprintf("%s ~ Virulence + (1|Isolate)", raw)), data = pa_cyto)
  p_lin = anova(lin, ddf = "Kenward-Roger")$`Pr(>F)`

  #estimated mean at each virulence score, and pairwise comparisons
  catm  = lmer(as.formula(sprintf("%s ~ factor(Virulence) + (1|Isolate)", raw)), data = pa_cyto)
  emm   = emmeans(catm, ~ Virulence)
  pairs_tab = as.data.frame(summary(pairs(emm, infer = c(TRUE, TRUE))))
  emm   = as.data.frame(emm)

  plot = ggplot(emm, aes(Virulence, emmean)) +
    geom_line() +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, linewidth = 0.8) +
    #isolate means
    geom_point(data = pa_cyto_iso, aes(x = Virulence, y = .data[[avg]]),
               inherit.aes = FALSE, size = 2,
               position = position_jitter(width = 0.1)) +
    #model means, drawn last so they sit on top
    geom_point(size = 5, shape = 21, fill = "black", colour = "white", stroke = 1.2) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = ANNOT_SZ,
             label = paste("Linear trend: p", fmt_p(p_lin))) +
    coord_cartesian(ylim = c(1, 4)) +
    labs(x = "Virulence Score", y = sprintf("Log10 %s (pg/mL)", CYTO$name[i])) +
    theme_fig()

  list(plot = plot, p_linear = p_lin, pairs = pairs_tab)
}

A = lapply(seq_len(nrow(CYTO)), trend_panel)

#linear-trend p-values and pairwise comparisons between virulence scores
tibble(cytokine = CYTO$name, p_linear = map_dbl(A, "p_linear"))
setNames(map(A, "pairs"), CYTO$name)


## ---- Panel B: Cox hazard ratios within dose quartiles -----------------

hr_panel = function(i) {
  avg = CYTO$avg[i]

  fit = coxph(as.formula(sprintf("Surv(Time, Dead) ~ %s * factor(cat_dose2)", avg)),
              data = mouse_mod)

  #overall cytokine x dose interaction
  p_int = car::Anova(fit, test.statistic = "Wald")[
    sprintf("%s:factor(cat_dose2)", avg), "Pr(>Chisq)"]

  #cytokine slope (log HR) within each quartile; Holm-adjusted test of HR = 1
  trend = emtrends(fit, ~ cat_dose2, var = avg, data = mouse_mod)
  tt    = as.data.frame(test(trend, adjust = "holm"))
  tr    = as.data.frame(trend)

  res = tibble(
    cat_dose = tr$cat_dose2,
    HR  = exp(tr[[paste0(avg, ".trend")]]),
    LCL = exp(tr$asymp.LCL),
    UCL = exp(tr$asymp.UCL),
    p   = tt$p.value
  ) %>%
    mutate(cat_dose = sub(",", "–", gsub("[][()]", "", cat_dose)),  # [5.1,6.6) -> 5.1–6.6
           cat_dose = factor(cat_dose, levels = cat_dose),
           label    = sprintf("%.2f (%.2f–%.2f) \np %s", HR, LCL, UCL, fmt_p(p)))

  plot = ggplot(res, aes(x = HR, y = cat_dose)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "gray90") +
    geom_errorbar(aes(xmin = LCL, xmax = UCL), width = 0.25, orientation = "y") +
    geom_point(size = 3) +
    #HR (95% CI) and p to the right of each bar; nudge_x is in HR units
    geom_text(aes(x = UCL, label = label), hjust = 0, nudge_x = 0.15, size = ANNOT_SZ) +
    scale_x_continuous(breaks = 0:5) +
    coord_cartesian(xlim = c(0, 5), clip = "off") +
    labs(x = sprintf("Hazard Ratio for %s", CYTO$name[i]),
         y = "Dose Quartile (log10 CFU)") +
    theme_fig()

  list(plot = plot, res = res, p_int = p_int, n = fit$n, nevent = fit$nevent,
       slope_pairs = pairs(trend, adjust = "holm"))
}

B = lapply(seq_len(nrow(CYTO)), hr_panel)

#hazard ratios, interaction p-values, model size, and slope comparisons
#between quartiles
setNames(map(B, "res"), CYTO$name)
tibble(cytokine = CYTO$name, p_interaction = map_dbl(B, "p_int"),
       n_mice = map_int(B, "n"), n_deaths = map_dbl(B, "nevent"))
setNames(map(B, "slope_pairs"), CYTO$name)


## ---- Panel C: Kaplan-Meier, low vs. high stimulators ------------------

#strains split at the median isolate level, among isolates in the mouse data
strain_lv = mouse_mod %>% distinct(ID, .keep_all = TRUE)

#survfit -> data frame, with an origin row and step-expanded CI bounds
km_tidy = function(d) {
  s   = summary(survfit(Surv(Time, Dead) ~ grp, data = d))
  out = data.frame(time = s$time, surv = s$surv, lo = s$lower, hi = s$upper,
                   grp = sub("grp=", "", as.character(s$strata)))
  org = data.frame(time = 0, surv = 1, lo = 1, hi = 1, grp = unique(out$grp))
  #carry each curve flat out to the end of follow-up
  end = out %>% group_by(grp) %>% slice_max(time, n = 1) %>%
    mutate(time = max(d$Time)) %>% ungroup() %>% filter(time > 0)
  out = bind_rows(org, out, end) %>% distinct() %>% arrange(grp, time)
  #duplicate each point at the next time so the ribbon steps with the curve
  out %>% group_by(grp) %>%
    reframe(time = c(time, lead(time)[-n()]),
            surv = c(surv, surv[-n()]),
            lo   = c(lo,   lo[-n()]),
            hi   = c(hi,   hi[-n()])) %>%
    arrange(grp, time, desc(surv))
}

#add the Low/High group for one cytokine
split_grp = function(d, v) {
  med = median(strain_lv[[v]])
  d %>% mutate(grp = factor(ifelse(.data[[v]] >= med, "High", "Low"),
                            levels = c("Low", "High")))
}

km_all = NULL; km_ann = NULL
for (i in seq_len(nrow(CYTO))) {
  for (qq in levels(mouse_mod$cat_dose2)) {
    d  = split_grp(filter(mouse_mod, cat_dose2 == qq), CYTO$avg[i])
    lr = survdiff(Surv(Time, Dead) ~ grp, data = d)
    km_all = bind_rows(km_all, km_tidy(d) %>% mutate(cyto = CYTO$expr[i], q = qq))
    km_ann = bind_rows(km_ann, data.frame(
      cyto = CYTO$expr[i], q = qq, p = 1 - pchisq(lr$chisq, 1),
      n = sprintf("Low n = %d, High n = %d", sum(d$grp == "Low"), sum(d$grp == "High"))))
  }
}

q_levels = pretty_bins(levels(mouse_mod$cat_dose2))
tidy_km  = function(df) df %>%
  mutate(cyto = factor(cyto, levels = CYTO$expr),
         q    = factor(pretty_bins(q), levels = q_levels))

#Holm adjustment across all 12 log-rank tests (3 cytokines x 4 quartiles)
km_ann = km_ann %>%
  mutate(p_adj = p.adjust(p, method = "holm"), lab = fmt_p_expr(p_adj)) %>%
  tidy_km()
km_all = tidy_km(km_all) %>% mutate(grp = factor(grp_lab[grp], levels = grp_lab))

#raw and adjusted p-values (supplementary table)
km_ann %>% select(cyto, q, p, p_adj) %>% arrange(cyto, q)

#shared KM layers; text size and colour differ between the two figures
km_layers = function(ann, txt_size, txt_col) list(
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.14, colour = NA),
  geom_line(linewidth = 0.7),
  geom_text(data = ann, aes(x = 0, y = 0.10, label = lab), parse = TRUE,
            inherit.aes = FALSE, hjust = 0, size = txt_size[1], colour = txt_col[1]),
  geom_text(data = ann, aes(x = 0, y = 0.02, label = n),
            inherit.aes = FALSE, hjust = 0, size = txt_size[2], colour = txt_col[2]),
  scale_colour_manual(values = setNames(c(col_low, col_high), grp_lab), name = NULL),
  scale_fill_manual(values = setNames(c(col_low, col_high), grp_lab), name = NULL),
  scale_x_continuous(breaks = c(0, 24, 48, 72, 96), limits = c(0, 96)),
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)),
  labs(x = "Time (hours)", y = "Survival")
)

#supplementary figure: every cytokine x dose quartile
km_plot = ggplot(km_all, aes(time, surv, colour = grp, fill = grp)) +
  km_layers(km_ann, txt_size = c(3.1, 2.7), txt_col = c("grey20", "grey40")) +
  facet_grid(cyto ~ q, labeller = labeller(cyto = label_parsed)) +
  theme_journal() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(colour = "grey10"))

ggsave(km_plot, file = "km_dose_quartiles.pdf", width = 13, height = 8)

#Figure 3 panel C: the focus quartile, one plot per cytokine
km_panel = function(i) {
  sel = function(df) filter(df, q == FOCUS_Q, cyto == CYTO$expr[i])
  ggplot(sel(km_all), aes(time, surv, colour = grp, fill = grp)) +
    km_layers(sel(km_ann), txt_size = c(ANNOT_SZ, ANNOT_SZ), txt_col = c("black", "black")) +
    labs(title = parse(text = CYTO$expr[i])) +
    theme_fig()
}

C = lapply(seq_len(nrow(CYTO)), km_panel)


## ---- Table: mortality by stimulation group within each dose bin -------
#DOSE_VAR = "cat_dose2" for the quartiles used in the figures, or
#"cat_dose" for the three nominal targets (6.5 / 7 / 7.75) used in Table 1.
DOSE_VAR = "cat_dose2"

km_table = NULL
for (i in seq_len(nrow(CYTO))) {
  for (qq in sort(unique(na.omit(mouse_mod[[DOSE_VAR]])))) {
    d = split_grp(filter(mouse_mod, .data[[DOSE_VAR]] == qq), CYTO$avg[i])
    if (n_distinct(d$grp) < 2) next
    lr = survdiff(Surv(Time, Dead) ~ grp, data = d)
    s  = d %>% group_by(grp) %>% summarise(n = n(), mort = mean(Dead), .groups = "drop")
    gv = function(k, col) s[[col]][s$grp == k]
    km_table = bind_rows(km_table, data.frame(
      Cytokine         = CYTO$name[i],
      Dose             = pretty_bins(as.character(qq)),
      Low_Stimulating  = sprintf("%.0f%% [%d]", 100 * gv("Low",  "mort"), gv("Low",  "n")),
      High_Stimulating = sprintf("%.0f%% [%d]", 100 * gv("High", "mort"), gv("High", "n")),
      LogRank_X2       = sprintf("%.2f", lr$chisq),
      P_value          = fmt_p_tab(1 - pchisq(lr$chisq, 1))))
  }
}

print(km_table, row.names = FALSE)
write.csv(km_table, "Table1_km_table.csv", row.names = FALSE, fileEncoding = "UTF-8")


## ---- Figure 3 assembly ------------------------------------------------

panels = c(map(A, "plot"), map(B, "plot"), C)
panels[[1]] = panels[[1]] + labs(tag = "A")
panels[[4]] = panels[[4]] + labs(tag = "B")
panels[[7]] = panels[[7]] + labs(tag = "C")

#the HR labels in row B sit to the right of each CI bar (clip = "off");
#extra right margin keeps them out of the next column. Panels stay equal.
panels[4:6] = lapply(panels[4:6], function(p) p + theme(plot.margin = margin(8, 70, 8, 8)))

#one flat 3 x 3 grid so every axis aligns and all panels are the same
#size; the KM legend is collected into one legend at the bottom
combined = wrap_plots(panels, ncol = 3, byrow = TRUE) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(combined, file = "Fig3_CombinedFig.pdf", width = 18, height = 16, device = cairo_pdf)
