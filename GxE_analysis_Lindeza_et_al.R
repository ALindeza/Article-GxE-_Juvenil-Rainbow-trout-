# =============================================================================
# A major life-history locus underlies genotype-by-environment variation
# in growth across water temperatures
#
# Lindeza et al.
# Analysis script — rainbow trout (Oncorhynchus mykiss) common-garden
# warming experiment, six6 genotype x temperature
#
# Data: Measurements_CorinneMasters_Final.xlsx  (n = 813 rows; 807 complete)
# =============================================================================


library(readxl)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(FactoMineR)
library(factoextra)
library(vegan)

data_path <- "Measurements_CorinneMasters_Final.xlsx"

raw <- read_excel(data_path)

str(raw)
colnames(raw)

# Six rows have no phenotype or genotype data and are dropped by every model - (n = 807 used throughout).
bad_rows <- which(!(is.finite(suppressWarnings(as.numeric(raw$Weight_Nov))) &
                      !is.na(raw$Treatment) & !is.na(raw$Six6) &
                      !is.na(raw$Family)   & !is.na(raw$Tank)))
bad_rows   # rows 808-813

#Treatment levels: C = control, T = warm (+2 degC).
# six6 heterozygotes are coded "H" in the raw file and relabelled "EL".
dat <- raw %>%
  mutate(
    Treatment = factor(Treatment, levels = c("C", "T")),
    Six6      = factor(ifelse(as.character(Six6) %in% c("H", "EL"), "EL",
                              as.character(Six6)),
                       levels = c("EE", "EL", "LL")),
    Sex       = factor(Sex),
    Family    = factor(Family),
    Tank      = factor(Tank)
  )

table(dat$Treatment, useNA = "ifany")
table(dat$Six6,      useNA = "ifany")

treat_cols <- c("C" = "#2C7FB8", "T" = "#D7301F")
geno_cols  <- c("EE" = "#08306B", "EL" = "#2171B5", "LL" = "#6BAED6")


# =============================================================================
# 1. TREATMENT x SEX MODELS  (Results 3.1)
#    Weight, length and condition at each of the three sampling points.
# =============================================================================

sex_models <- list(
  Weight_Nov    = Weight_Nov    ~ Treatment * Sex + (1 | Family) + (1 | Tank),
  Weight_Mar    = Weight_Mar    ~ Treatment * Sex + (1 | Family) + (1 | Tank),
  Weight_May    = Weight_May    ~ Treatment * Sex + (1 | Family) + (1 | Tank),
  Length_Nov    = Length_Nov    ~ Treatment * Sex + (1 | Family) + (1 | Tank),
  Length_Mar    = Length_Mar    ~ Treatment * Sex + (1 | Family) + (1 | Tank),
  Length_May    = Length_May    ~ Treatment * Sex + (1 | Family) + (1 | Tank),
  Condition_Nov = Condition_Nov ~ Treatment * Sex + (1 | Family) + (1 | Tank),
  Condition_Mar = Condition_Mar ~ Treatment * Sex + (1 | Family) + (1 | Tank),
  Condition_May = Condition_May ~ Treatment * Sex + (1 | Family) + (1 | Tank)
)

fits_sex <- lapply(sex_models, function(f) lmer(f, data = dat))

lapply(fits_sex, anova)     # Table S2
lapply(fits_sex, summary)

# Family as random vs fixed (May weight), fitted by ML for a valid AIC comparison
AIC(
  lmer(Weight_May ~ Treatment * Sex + (1 | Family) + (1 | Tank),
       data = dat, REML = FALSE),
  lmer(Weight_May ~ Treatment * Sex + Family + (1 | Tank),
       data = dat, REML = FALSE)
)


# =============================================================================
# 2. TREATMENT x six6 MODELS  (Results 3.2)
# =============================================================================

# Nov and Mar models, and all condition models, give singular fits with Tank
# included (Tank variance estimated at 0); those are refitted with Family only.
geno_models <- list(
  Weight_Nov    = Weight_Nov    ~ Treatment * Six6 + (1 | Family) + (1 | Tank),
  Weight_Mar    = Weight_Mar    ~ Treatment * Six6 + (1 | Family) + (1 | Tank),
  Weight_May    = Weight_May    ~ Treatment * Six6 + (1 | Family),
  Length_Nov    = Length_Nov    ~ Treatment * Six6 + (1 | Family) + (1 | Tank),
  Length_Mar    = Length_Mar    ~ Treatment * Six6 + (1 | Family) + (1 | Tank),
  Length_May    = Length_May    ~ Treatment * Six6 + (1 | Family),
  Condition_Nov = Condition_Nov ~ Treatment * Six6 + (1 | Family),
  Condition_Mar = Condition_Mar ~ Treatment * Six6 + (1 | Family),
  Condition_May = Condition_May ~ Treatment * Six6 + (1 | Family)
)

fits_geno <- lapply(geno_models, function(f) lmer(f, data = dat))

lapply(fits_geno, anova)    # Table S3
lapply(fits_geno, summary)

m_weight_may <- fits_geno$Weight_May
m_length_may <- fits_geno$Length_May

# --- Estimated marginal means and contrasts, May body weight -----------------
emm_w_by_geno  <- emmeans(m_weight_may, ~ Treatment | Six6)   # reaction norms
emm_w_by_treat <- emmeans(m_weight_may, ~ Six6 | Treatment)   # genotype effects

as.data.frame(emm_w_by_geno)
summary(pairs(emm_w_by_geno),  infer = TRUE)   # warm - control within genotype
summary(pairs(emm_w_by_treat), infer = TRUE)   # genotype contrasts within env.
summary(contrast(emm_w_by_geno, interaction = "pairwise"), infer = TRUE)

# Same contrasts for May body length
summary(pairs(emmeans(m_length_may, ~ Treatment | Six6)), infer = TRUE)
summary(pairs(emmeans(m_length_may, ~ Six6 | Treatment)), infer = TRUE)

# --- Size-corrected weight model (length as covariate) ----------------------  # Table S4
m_weight_may_adj <- lmer(Weight_May ~ Treatment * Six6 + Length_May + (1 | Family),
                         data = dat)
anova(m_weight_may_adj)
summary(m_weight_may_adj)

emm_adj <- emmeans(m_weight_may_adj, ~ Treatment | Six6)
summary(pairs(emm_adj), infer = TRUE)


# =============================================================================
# 3. BODY COMPOSITION — DXA TRAITS  (Results 3.3)
# =============================================================================

dxa_traits <- c("Fat", "Fat_perc", "Lean", "BMD", "BMC")

fits_dxa <- lapply(setNames(dxa_traits, dxa_traits), function(tr) {
  lmer(as.formula(paste(tr, "~ Treatment * Six6 + (1 | Family)")), data = dat)
})

lapply(fits_dxa, anova)     # Table S5
lapply(fits_dxa, summary)

# Genotype contrasts within each thermal environment                            # Table S6
lapply(fits_dxa, function(m) emmeans(m, pairwise ~ Six6 | Treatment))

# --- Size-adjusted composition models (log body weight as covariate) --------  # Table S7
dat$logW <- as.numeric(scale(log(dat$Weight_May), center = TRUE, scale = FALSE))

fits_dxa_adj <- lapply(setNames(c("Lean", "Fat", "BMC"), c("Lean", "Fat", "BMC")),
                       function(tr) {
                         lmer(as.formula(paste(tr,
                              "~ Treatment * Six6 + logW + (1 | Family) + (1 | Tank)")),
                              data = dat)
                       })

lapply(fits_dxa_adj, function(m) {
  pairs(emmeans(m, ~ Treatment * Six6, cov.reduce = list(logW = mean)))
})

# =============================================================================
# 4. MULTIVARIATE BODY COMPOSITION — PCA AND PERMANOVA  (Results 3.3)
# =============================================================================

# Residualise each DXA trait on final body weight, with family as a random
# effect, to remove overall size.
resid_trait <- function(trait, size_var = "Weight_May") {
  f <- as.formula(paste(trait, "~", size_var, "+ (1 | Family)"))
  resid(lmer(f, data = dat, na.action = na.exclude))
}

pca_dat <- dat %>%
  mutate(
    Lean_r      = resid_trait("Lean"),
    Fat_r       = resid_trait("Fat"),
    Fat_perc_r  = resid_trait("Fat_perc"),
    BMD_r       = resid_trait("BMD"),
    BMC_r       = resid_trait("BMC"),
    Bone_area_r = resid_trait("Bone_area")
  ) %>%
  select(Lean_r, Fat_r, Fat_perc_r, BMD_r, BMC_r, Bone_area_r,
         Treatment, Family) %>%
  na.omit()

pca_resid <- PCA(select(pca_dat, ends_with("_r")),
                 scale.unit = TRUE, graph = FALSE)

pca_resid$eig          # Table S8a
pca_resid$var$cor      # Table S8b
pca_resid$var$contrib

fviz_pca_var(pca_resid, col.var = "contrib",
             gradient.cols = c("red4", "sienna1"), repel = TRUE) +
  ggtitle("PCA variable loadings - weight-corrected residuals") +
  theme_minimal()

fviz_pca_biplot(pca_resid,
                col.ind = pca_dat$Treatment,
                palette = c("C" = "steelblue", "T" = "firebrick"),
                col.var = "black", addEllipses = TRUE,
                ellipse.type = "norm", legend.title = "Treatment",
                repel = TRUE, label = "var") +
  labs(x = paste0("PC1 (", round(pca_resid$eig[1, 2], 1), "% variance)"),
       y = paste0("PC2 (", round(pca_resid$eig[2, 2], 1), "% variance)")) +
  theme_minimal()

# PERMANOVA on the first two PC scores, terms tested marginally                 # Table S9
scores_df <- as.data.frame(pca_resid$ind$coord[, 1:2])
scores_df$Treatment <- pca_dat$Treatment
scores_df$Family    <- factor(pca_dat$Family)

adonis2(scores_df[, 1:2] ~ Family + Treatment,
        data = scores_df, permutations = 999,
        method = "euclidean", by = "margin")


# =============================================================================
# 5. PROPORTION OF FAMILY VARIANCE (PFV)  (Results 3.4)
#
# Family-based variance partitioning, estimated separately per environment.
# With only 4 families this captures additive + dominance + common-environment
# effects and is interpreted as family-level resemblance, not narrow-sense h2.
# =============================================================================

fit_pfv <- function(df) {
  mod <- lmer(Weight_May ~ Sex + Six6 + (1 | Family) + (1 | Tank),
              data = df, REML = TRUE)
  vc <- as.data.frame(VarCorr(mod))
  V_family <- vc$vcov[vc$grp == "Family"]
  V_tank   <- vc$vcov[vc$grp == "Tank"]
  V_resid  <- vc$vcov[vc$grp == "Residual"]
  list(model = mod, V_family = V_family, V_tank = V_tank, V_resid = V_resid,
       V_total = V_family + V_tank + V_resid,
       pfv = V_family / (V_family + V_tank + V_resid))
}

# Parametric bootstrap CIs.
# use.u = TRUE conditions on the observed random effects; the marginal
# bootstrap produces degenerate CIs with only 4 families.
boot_pfv <- function(fitted, nsim = 999, seed = 123) {
  set.seed(seed)
  stat_fun <- function(m) {
    vc <- as.data.frame(VarCorr(m))
    vc$vcov[vc$grp == "Family"] /
      (vc$vcov[vc$grp == "Family"] + vc$vcov[vc$grp == "Tank"] +
         vc$vcov[vc$grp == "Residual"])
  }
  bb <- bootMer(fitted$model, FUN = stat_fun, nsim = nsim,
                use.u = TRUE, type = "parametric", parallel = "no")
  list(boot = bb, ci = quantile(bb$t, c(0.025, 0.975), na.rm = TRUE))
}

# One-sided LRT (variance components are bounded at zero)
lrt_family <- function(fitted) {
  full   <- fitted$model
  no_fam <- update(full, . ~ . - (1 | Family))
  lrt    <- anova(no_fam, full, refit = FALSE)
  list(lrt = lrt, p_one_sided = lrt$`Pr(>Chisq)`[2] / 2)
}

dat_w <- filter(dat, !is.na(Weight_May))
res_C <- fit_pfv(filter(dat_w, Treatment == "C"))
res_T <- fit_pfv(filter(dat_w, Treatment == "T"))

lrt_C <- lrt_family(res_C); lrt_C$lrt; lrt_C$p_one_sided
lrt_T <- lrt_family(res_T); lrt_T$lrt; lrt_T$p_one_sided

ci_C <- boot_pfv(res_C)
ci_T <- boot_pfv(res_T)

# Difference between environments                                               # Table S10
delta_pfv <- ci_C$boot$t - ci_T$boot$t
round(quantile(delta_pfv, c(0.025, 0.975), na.rm = TRUE), 3)

pfv_summary <- tibble(
  Treatment = c("Control", "Warm"),
  V_family  = c(res_C$V_family, res_T$V_family),
  V_tank    = c(res_C$V_tank,   res_T$V_tank),
  V_resid   = c(res_C$V_resid,  res_T$V_resid),
  V_total   = c(res_C$V_total,  res_T$V_total),
  PFV       = c(res_C$pfv,      res_T$pfv),
  CI_lower  = c(ci_C$ci[1],     ci_T$ci[1]),
  CI_upper  = c(ci_C$ci[2],     ci_T$ci[2]),
  LRT_p     = c(lrt_C$p_one_sided, lrt_T$p_one_sided)
)
pfv_summary

summary(res_C$model)
summary(res_T$model)


# =============================================================================
# 6. FIGURES
# =============================================================================

# ---- 6.1 Trait distributions by treatment (one panel per timepoint) ---------
# Helper replacing twelve near-identical ggplot blocks.
violin_by <- function(x, y, ylab, title, data = dat) {
  d <- data %>% filter(!is.na(.data[[y]]), !is.na(.data[[x]]))
  ggplot(d, aes(x = .data[[x]], y = .data[[y]], fill = Treatment)) +
    geom_violin(alpha = 0.3, width = 0.9, trim = TRUE) +
    geom_boxplot(width = 0.25, position = position_dodge(width = 0.9),
                 outlier.shape = NA, alpha = 0.6) +
    geom_jitter(aes(color = Treatment),
                position = position_jitterdodge(jitter.width = 0.2,
                                                dodge.width = 0.8),
                size = 1.6, alpha = 0.6) +
    scale_fill_manual(values = treat_cols) +
    scale_color_manual(values = treat_cols) +
    labs(x = if (x == "Six6") "six6 genotype" else "Treatment",
         y = ylab, title = title) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0.5))
}

# By treatment only
violin_by("Treatment", "Weight_Nov", "Weight (g)",  "Weight in November")
violin_by("Treatment", "Weight_Mar", "Weight (g)",  "Weight in March")
violin_by("Treatment", "Weight_May", "Weight (g)",  "Weight in May")
violin_by("Treatment", "Length_Nov", "Length (cm)", "Length in November")
violin_by("Treatment", "Length_Mar", "Length (cm)", "Length in March")
violin_by("Treatment", "Length_May", "Length (cm)", "Length in May")

# By six6 genotype and treatment
violin_by("Six6", "Weight_Nov", "Weight (g)",  "Weight in November by six6")
violin_by("Six6", "Weight_Mar", "Weight (g)",  "Weight in March by six6")
violin_by("Six6", "Weight_May", "Weight (g)",  "Weight in May by six6")
violin_by("Six6", "Length_Nov", "Length (cm)", "Length in November by six6")
violin_by("Six6", "Length_Mar", "Length (cm)", "Length in March by six6")
violin_by("Six6", "Length_May", "Length (cm)", "Length in May by six6")

# DXA traits by six6 genotype and treatment
violin_by("Six6", "Fat",      "Fat (g)",  "Fat in May by six6")
violin_by("Six6", "Fat_perc", "Fat (%)",  "Fat (%) in May by six6")
violin_by("Six6", "Lean",     "Lean (g)", "Lean mass in May by six6")
violin_by("Six6", "BMC",      "BMC",      "BMC in May by six6")
violin_by("Six6", "BMD",      "BMD",      "BMD in May by six6")

# ---- 6.2 Density distributions by six6 genotype ----------------------------
density_by_geno <- function(y, xlab, title) {
  ggplot(dat, aes(x = .data[[y]], fill = Treatment)) +
    geom_density(alpha = 0.4) +
    facet_wrap(~ Six6) +
    scale_fill_manual(values = treat_cols) +
    labs(title = title, x = xlab, y = "Density") +
    theme_bw()
}

density_by_geno("Weight_Nov", "Weight (Nov, g)", "November weight by six6")
density_by_geno("Weight_Mar", "Weight (Mar, g)", "March weight by six6")
density_by_geno("Weight_May", "Weight (May, g)", "May weight by six6")

# ---- 6.3 Growth trajectories by genotype and environment -------------------
traj_data <- function(prefix, value_name) {
  dat %>%
    select(Treatment, Six6, all_of(paste0(prefix, c("_Nov", "_Mar", "_May")))) %>%
    pivot_longer(starts_with(prefix),
                 names_to = "TimePoint", values_to = value_name) %>%
    mutate(TimePoint = factor(sub(".*_", "", TimePoint),
                              levels = c("Nov", "Mar", "May"),
                              labels = c("November", "March", "May"))) %>%
    filter(!is.na(.data[[value_name]]), !is.na(Treatment), !is.na(Six6)) %>%
    group_by(Treatment, Six6, TimePoint) %>%
    summarise(mean_val = mean(.data[[value_name]]),
              se_val   = sd(.data[[value_name]]) / sqrt(n()),
              n = n(), .groups = "drop")
}

traj_plot <- function(sum_df, ylab) {
  ggplot(sum_df, aes(x = TimePoint, y = mean_val, group = Six6, color = Six6)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                  width = 0.08, linewidth = 0.8) +
    facet_wrap(~ Treatment, nrow = 1) +
    scale_color_manual(values = geno_cols) +
    labs(x = "Time point", y = ylab, color = "Genotype") +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold", size = 13),
          axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "top",
          axis.title = element_text(face = "bold"))
}

traj_plot(traj_data("Weight", "Weight"), "Mean body weight (g)")
traj_plot(traj_data("Length", "Length"), "Mean body length (cm)")

# ---- 6.4 Reaction norms (model-based) --------------------------------------
# NOTE: these use the same model specification as Section 2 (Family + Tank
# where estimable), so that plotted slopes match the reported contrasts.
rxn_plot <- function(model, title) {
  emmd <- as.data.frame(emmeans(model, ~ Treatment | Six6))
  ggplot(emmd, aes(x = Treatment, y = emmean, group = Six6, color = Six6)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.4) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                  width = 0.05, linewidth = 0.6) +
    scale_color_manual(values = geno_cols, name = "six6") +
    labs(title = title, x = "Environment (C -> T)",
         y = "Adjusted mean weight (g)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0.5))
}

rxn_plot(fits_geno$Weight_Nov, "Reaction norms - November")
rxn_plot(fits_geno$Weight_Mar, "Reaction norms - March")
rxn_plot(fits_geno$Weight_May, "Reaction norms - May")

# Slopes (warm - control) per genotype
lapply(fits_geno[c("Weight_Nov", "Weight_Mar", "Weight_May")], function(m) {
  summary(contrast(emmeans(m, ~ Treatment | Six6),
                   method = "revpairwise", by = "Six6"), infer = TRUE)
})

# ---- 6.5 Variance partitioning and PFV -------------------------------------
var_table <- tibble(
  Environment = rep(c("Control", "Warm"), each = 3),
  Component   = factor(rep(c("Family", "Tank", "Residual"), 2),
                       levels = c("Residual", "Tank", "Family")),
  Variance    = c(res_C$V_family, res_C$V_tank, res_C$V_resid,
                  res_T$V_family, res_T$V_tank, res_T$V_resid)
)

ggplot(var_table, aes(x = Environment, y = Variance, fill = Component)) +
  geom_col(position = "stack", width = 0.5) +
  scale_fill_manual(values = c("Family" = "#2c7bb6", "Tank" = "#abd9e9",
                               "Residual" = "#d9d9d9")) +
  labs(title = "Variance partitioning of body weight by thermal environment",
       x = NULL, y = "Variance estimate") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(pfv_summary, aes(x = Treatment, y = PFV, fill = Treatment)) +
  geom_col(width = 0.45) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                width = 0.12, linewidth = 0.8) +
  scale_fill_manual(values = c("Control" = "steelblue", "Warm" = "firebrick")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Proportion of family variance in body weight",
       subtitle = "Parametric bootstrap 95% CI",
       x = NULL, y = "PFV") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 10, colour = "grey40"))

