# ============================================================
# CNA × Treatment Interaction Analysis — HNSCC
# MSc Cancer and Clinical Oncology Dissertation
# Author: Mrinalini Sarkar
# Date: June 2026
#
# Key findings:
# - CCND1 amplification: HR=2.89 on CDK4/6i (p=0.014, BH=0.049)
# - CDKN2A copy number loss: HR=0.234 on CDK4/6i (p=0.002, BH=0.011)
# - CDKN2A × Immunotherapy: HR=6.28 (p=0.002, BH=0.011) — reversal
# ============================================================

# ── 0. LIBRARIES ─────────────────────────────────────────────
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(survival)
library(broom)
library(ggplot2)
library(writexl)
library(survminer)

# ── 1. LOAD DATA ─────────────────────────────────────────────
master_raw <- read_excel(
  "CNA_Clean_Master (2) copy.xlsx",
  sheet     = "Master Matrix",
  col_names = FALSE
)

# Row 2 contains actual column headers
# Row 1 contains section banners — skip both
colnames(master_raw) <- as.character(master_raw[2, ])
master <- master_raw[-c(1, 2), ]

# ── 2. RENAME AND CONVERT COLUMNS ────────────────────────────
master <- master %>%
  rename(
    patient_id       = `Patient ID`,
    study            = `Study`,
    cancer_type      = `Cancer\r\nType`,
    treatment_class  = `Treatment Class`,
    drug             = `Drug`,
    setting          = `Setting`,
    hpv              = `HPV`,
    bor              = `BOR`,
    clinical_benefit = `Clinical\r\nBenefit`,
    pfs_months       = `PFS\r\n(mo)`,
    os_months        = `OS\r\n(mo)`,
    cna_notes        = `CNA Notes`
  ) %>%
  mutate(
    pfs_months = as.numeric(pfs_months),
    os_months  = as.numeric(os_months)
  )

cat("Rows loaded:", nrow(master), "\n")

# ── 3. CNA BINARY VARIABLES ──────────────────────────────────
# IMPORTANT: Mutations (Mut) are NOT CNAs
# Only Amplification (Amp) and Deletion (Del) count
# TP53 (Del=2) and TERT (Amp=1) excluded — too rare for CNA analysis

master <- master %>%
  mutate(
    CCND1_cna  = ifelse(CCND1  == "Amp", 1, 0),
    CDKN2A_cna = ifelse(CDKN2A %in% c("Amp", "Del"), 1, 0),
    EGFR_cna   = ifelse(EGFR   %in% c("Amp", "Del"), 1, 0),
    PIK3CA_cna = ifelse(PIK3CA == "Amp", 1, 0),
    pfs_event  = ifelse(is.na(pfs_months), 0, 1)
  )

# Verify — expected: CCND1=61, CDKN2A=23, EGFR=39, PIK3CA=23
cat("\n=== CNA COUNTS (should be 61/23/39/23) ===\n")
cat("CCND1_cna:",  sum(master$CCND1_cna,  na.rm=TRUE), "\n")
cat("CDKN2A_cna:", sum(master$CDKN2A_cna, na.rm=TRUE), "\n")
cat("EGFR_cna:",   sum(master$EGFR_cna,   na.rm=TRUE), "\n")
cat("PIK3CA_cna:", sum(master$PIK3CA_cna, na.rm=TRUE), "\n")

# ── 3b. CLINICAL BENEFIT VARIABLE ────────────────────────────
master <- master %>%
  mutate(
    benefit = case_when(
      clinical_benefit == "CB"  ~ "Benefit",
      clinical_benefit == "NCB" ~ "No Benefit",
      TRUE ~ NA_character_
    )
  )

# ── 4. TREATMENT GROUP AND HPV HARMONISATION ─────────────────
master <- master %>%
  mutate(
    Treatment_Group = case_when(
      grepl("CDK4/6",        treatment_class) ~ "CDK4/6",
      grepl("EGFR",          treatment_class) ~ "EGFR",
      grepl("Immunotherapy", treatment_class) ~ "Immunotherapy",
      grepl("PI3K",          treatment_class) ~ "PI3K",
      TRUE ~ "Other"
    ),
    # Harmonise HPV from 8+ coding schemes to 3 categories
    HPV_clean = case_when(
      hpv %in% c("Positive", "16", "33") ~ "Positive",
      hpv %in% c("Negative", "N", "WT")  ~ "Negative",
      TRUE ~ "Unknown"
    )
  )

cat("\n=== TREATMENT GROUP COUNTS ===\n")
print(table(master$Treatment_Group))

# ── 5. COHORT OVERVIEW ───────────────────────────────────────
cohort_overview <- master %>%
  group_by(Treatment_Group) %>%
  summarise(
    N              = n(),
    BOR_available  = sum(bor %in% c("CR","PR","SD","PD","NE")),
    CB_evaluable   = sum(bor %in% c("CR","PR","SD","PD")),
    PFS_available  = sum(!is.na(pfs_months)),
    OS_available   = sum(!is.na(os_months)),
    CCND1_CNA_n    = sum(CCND1_cna,  na.rm=TRUE),
    CDKN2A_CNA_n   = sum(CDKN2A_cna, na.rm=TRUE),
    EGFR_CNA_n     = sum(EGFR_cna,   na.rm=TRUE),
    PIK3CA_CNA_n   = sum(PIK3CA_cna, na.rm=TRUE),
    .groups = "drop"
  )

cat("\n=== COHORT OVERVIEW ===\n")
print(cohort_overview)

# ── 6. UNIVARIATE COX MODELS ─────────────────────────────────
genes_cna <- c("CCND1_cna", "CDKN2A_cna", "EGFR_cna", "PIK3CA_cna")

uni_results <- lapply(genes_cna, function(g) {
  fit <- coxph(
    as.formula(paste0("Surv(pfs_months) ~ ", g)),
    data = master
  )
  tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>%
    mutate(gene = g)
}) %>%
  bind_rows() %>%
  select(gene, estimate, conf.low, conf.high, p.value) %>%
  rename(HR=estimate, CI_low=conf.low, CI_high=conf.high, p=p.value)

cat("\n=== UNIVARIATE COX (expect all null — confirms treatment context needed) ===\n")
print(uni_results)

# ── 7. MULTIVARIABLE INTERACTION MODELS ──────────────────────
# CDK4/6 = reference group
# HPV_clean = covariate
# Exclude Other group (no outcome data)

master_model <- master %>%
  filter(Treatment_Group != "Other",
         !is.na(pfs_months)) %>%
  mutate(
    Treatment_Group = relevel(factor(Treatment_Group), ref = "CDK4/6")
  )

cat("\n=== FITTING INTERACTION MODELS ===\n")

cox_ccnd1 <- coxph(
  Surv(pfs_months) ~ CCND1_cna * Treatment_Group + HPV_clean,
  data = master_model
)

cox_cdkn2a <- coxph(
  Surv(pfs_months) ~ CDKN2A_cna * Treatment_Group + HPV_clean,
  data = master_model
)

cox_egfr <- coxph(
  Surv(pfs_months) ~ EGFR_cna * Treatment_Group + HPV_clean,
  data = master_model
)

ccnd1_results  <- tidy(cox_ccnd1,  exponentiate=TRUE, conf.int=TRUE)
cdkn2a_results <- tidy(cox_cdkn2a, exponentiate=TRUE, conf.int=TRUE)
egfr_results   <- tidy(cox_egfr,   exponentiate=TRUE, conf.int=TRUE)

cat("\n--- CCND1 model ---\n")
print(ccnd1_results %>%
        filter(grepl("CCND1", term)) %>%
        select(term, estimate, conf.low, conf.high, p.value))

cat("\n--- CDKN2A model ---\n")
print(cdkn2a_results %>%
        filter(grepl("CDKN2A", term)) %>%
        select(term, estimate, conf.low, conf.high, p.value))

cat("\n--- EGFR model ---\n")
print(egfr_results %>%
        filter(grepl("EGFR_cna", term)) %>%
        select(term, estimate, conf.low, conf.high, p.value))

# ── 8. COMBINED HR PER TREATMENT GROUP (exact term matching) ─────────────────

combined_hr <- tibble(
  Gene = c(
    "CCND1 Amplification",
    "CCND1 Amplification",
    "CCND1 Amplification",
    "CDKN2A Copy Number Loss",
    "CDKN2A Copy Number Loss"
  ),
  Treatment = c(
    "CDK4/6 Inhibitor (reference)",
    "EGFR-targeted",
    "Immunotherapy",
    "CDK4/6 Inhibitor (reference)",
    "Immunotherapy"
  ),
  # Main effect HR (direct from CDK4/6i reference row)
  Main_effect_HR = c(
    get_est(ccnd1_results,  "CCND1_cna"),
    get_est(ccnd1_results,  "CCND1_cna"),
    get_est(ccnd1_results,  "CCND1_cna"),
    get_est(cdkn2a_results, "CDKN2A_cna"),
    get_est(cdkn2a_results, "CDKN2A_cna")
  ),
  # Interaction HR (1.0 for reference group — no interaction term)
  Interaction_HR = c(
    1.0,
    get_est(ccnd1_results,  "CCND1_cna:Treatment_GroupEGFR"),
    get_est(ccnd1_results,  "CCND1_cna:Treatment_GroupImmunotherapy"),
    1.0,
    get_est(cdkn2a_results, "CDKN2A_cna:Treatment_GroupImmunotherapy")
  ),
  # Interaction p-value
  Interaction_p = c(
    get_p(ccnd1_results,  "CCND1_cna"),
    get_p(ccnd1_results,  "CCND1_cna:Treatment_GroupEGFR"),
    get_p(ccnd1_results,  "CCND1_cna:Treatment_GroupImmunotherapy"),
    get_p(cdkn2a_results, "CDKN2A_cna"),
    get_p(cdkn2a_results, "CDKN2A_cna:Treatment_GroupImmunotherapy")
  )
) %>%
  mutate(
    # Combined HR = main effect × interaction term
    Combined_HR = round(Main_effect_HR * Interaction_HR, 3),
    Interpretation = c(
      "Harmful — 3x hazard of progression on CDK4/6i",
      "Effect attenuated — near neutral on EGFR-targeted",
      "Effect attenuated — not significant on immunotherapy",
      "Protective — 77% lower hazard on CDK4/6i",
      "Reversed — harmful on immunotherapy (HPV-adjusted)"
    )
  ) %>%
  select(Gene, Treatment, Combined_HR, Interaction_p, Interpretation)

cat("\n=== COMBINED HR BY TREATMENT GROUP ===\n")
print(combined_hr)

# ── 9. BH MULTIPLE TESTING CORRECTION ───────────────────────
all_p <- tibble(
  test = c(
    "Univariate: CCND1_cna",
    "Univariate: CDKN2A_cna",
    "Univariate: EGFR_cna",
    "Univariate: PIK3CA_cna",
    "CCND1_cna — CDK4/6i (main effect)",
    "CCND1_cna × EGFR-targeted",
    "CCND1_cna × Immunotherapy",
    "CDKN2A_cna — CDK4/6i (main effect)",
    "CDKN2A_cna × Immunotherapy",
    "EGFR_cna — CDK4/6i (main effect)",
    "EGFR_cna × EGFR-targeted"
  ),
  p_raw = c(
    uni_results$p[1],
    uni_results$p[2],
    uni_results$p[3],
    uni_results$p[4],
    get_p(ccnd1_results,  "CCND1_cna$"),
    get_p(ccnd1_results,  "EGFR"),
    get_p(ccnd1_results,  "Immunotherapy"),
    get_p(cdkn2a_results, "CDKN2A_cna$"),
    get_p(cdkn2a_results, "Immunotherapy"),
    get_p(egfr_results,   "EGFR_cna$"),
    get_p(egfr_results,   "EGFR_cna:Treatment")
  )
) %>%
  filter(!is.na(p_raw)) %>%
  mutate(
    p_BH     = p.adjust(p_raw, method = "BH"),
    sig_raw  = p_raw < 0.05,
    sig_BH   = p_BH  < 0.05
  ) %>%
  arrange(p_raw)

cat("\n=== BH CORRECTION RESULTS ===\n")
print(all_p)
cat("\nTests surviving BH correction:\n")
print(all_p %>% filter(sig_BH))

# ── 10. DESCRIPTIVE SUMMARY TABLE ────────────────────────────
cna_summary <- master %>%
  filter(Treatment_Group != "Other") %>%
  pivot_longer(
    cols      = c(CCND1_cna, CDKN2A_cna, EGFR_cna),
    names_to  = "Gene",
    values_to = "CNA_status"
  ) %>%
  mutate(
    Gene = recode(Gene,
                  "CCND1_cna"  = "CCND1 Amplification",
                  "CDKN2A_cna" = "CDKN2A Copy Number Loss",
                  "EGFR_cna"   = "EGFR Copy Number Alteration"
    ),
    CNA_label = ifelse(CNA_status == 1, "CNA+", "CNA-")
  ) %>%
  group_by(Gene, Treatment_Group, CNA_label) %>%
  summarise(
    n          = n(),
    CB_n       = sum(benefit == "Benefit", na.rm=TRUE),
    NCB_n      = sum(benefit == "No Benefit", na.rm=TRUE),
    CB_pct     = round(100 * CB_n / (CB_n + NCB_n), 1),
    median_pfs = round(median(pfs_months, na.rm=TRUE), 1),
    n_pfs      = sum(!is.na(pfs_months)),
    .groups    = "drop"
  ) %>%
  filter(n >= 3) %>%
  arrange(Gene, Treatment_Group, CNA_label)

cat("\n=== DESCRIPTIVE SUMMARY ===\n")
print(cna_summary, n=50)

# ── 11. FIGURE A: FOREST PLOT ────────────────────────────────
forest_data <- tibble(
  term  = c(
    "CCND1 Amp\n(CDK4/6i)",
    "CCND1 Amp\n× EGFR-targeted",
    "CCND1 Amp\n× Immunotherapy",
    "CDKN2A CNL\n(CDK4/6i)",
    "CDKN2A CNL\n× Immunotherapy",
    "EGFR CNA\n(CDK4/6i)"
  ),
  HR    = c(2.89, 0.336, 0.491, 0.234, 6.28, 0.586),
  lower = c(1.25, 0.118, 0.182, 0.0948, 1.96, 0.326),
  upper = c(6.70, 0.956, 1.32,  0.577,  20.1, 1.05),
  gene  = c("CCND1","CCND1","CCND1","CDKN2A","CDKN2A","EGFR"),
  p_BH  = c(0.0494, 0.113, 0.293, 0.0108, 0.0108, 0.164),
  sig   = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE)
)

ggplot(forest_data,
       aes(y = reorder(term, HR), x = HR,
           colour = gene, shape = sig)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), size = 0.7) +
  geom_vline(xintercept = 1, linetype="dashed", colour="grey40") +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 20)) +
  scale_colour_manual(values = c(
    "CCND1"  = "#2E75B6",
    "CDKN2A" = "#C0504D",
    "EGFR"   = "#70AD47"
  )) +
  scale_shape_manual(
    values = c("TRUE"=16, "FALSE"=1),
    labels = c("TRUE"="BH p<0.05", "FALSE"="BH p≥0.05")
  ) +
  labs(
    title    = "CNA × Treatment Interaction Effects on PFS in HNSCC",
    subtitle = "Filled = survives BH correction | Open = does not\nReference: CDK4/6 inhibitor, CNA wildtype",
    x        = "Hazard Ratio (log scale)",
    y        = "",
    colour   = "Gene",
    shape    = "BH significance"
  ) +
  theme_minimal(base_size=12) +
  theme(legend.position="bottom")

ggsave("Figure_A_Forest_CNA_interactions.png",
       width=9, height=6, dpi=300)
cat("Saved: Figure_A_Forest_CNA_interactions.png\n")

# ── 12. FIGURE B: CB RATE BY TREATMENT ───────────────────────
cb_plot <- cna_summary %>%
  filter(Gene != "EGFR Copy Number Alteration",
         !is.na(CB_pct),
         (CB_n + NCB_n) >= 3)

ggplot(cb_plot,
       aes(x = Treatment_Group, y = CB_pct, fill = CNA_label)) +
  geom_col(position="dodge", width=0.6) +
  geom_text(aes(label = paste0(CB_pct, "%\n(n=", n, ")")),
            position = position_dodge(0.6),
            vjust = -0.3, size = 2.8) +
  facet_wrap(~ Gene, ncol=1, scales="free_x") +
  scale_fill_manual(values=c("CNA+"="#E8836B", "CNA-"="#7FB37F")) +
  scale_y_continuous(limits=c(0,120),
                     labels=scales::percent_format(scale=1)) +
  labs(
    title    = "Clinical Benefit Rate by CNA Status and Treatment Class",
    subtitle = "CB = CR + PR + SD | Groups with n<3 excluded",
    x        = "Treatment Class",
    y        = "Clinical Benefit Rate (%)",
    fill     = "CNA Status"
  ) +
  theme_minimal(base_size=12) +
  theme(
    axis.text.x     = element_text(angle=40, hjust=1),
    strip.text      = element_text(face="bold"),
    legend.position = "top"
  )

ggsave("Figure_B_CB_rate_by_treatment.png",
       width=10, height=9, dpi=300)
cat("Saved: Figure_B_CB_rate_by_treatment.png\n")

# ── 13. FIGURE C: MEDIAN PFS BY TREATMENT ────────────────────
pfs_plot <- cna_summary %>%
  filter(Gene != "EGFR Copy Number Alteration",
         !is.na(median_pfs),
         n_pfs >= 3)

ggplot(pfs_plot,
       aes(x = Treatment_Group, y = median_pfs, fill = CNA_label)) +
  geom_col(position="dodge", width=0.6) +
  geom_text(aes(label = paste0(median_pfs, "m\n(n=", n_pfs, ")")),
            position = position_dodge(0.6),
            vjust = -0.3, size = 2.8) +
  facet_wrap(~ Gene, ncol=1, scales="free_x") +
  scale_fill_manual(values=c("CNA+"="#E8836B", "CNA-"="#7FB37F")) +
  labs(
    title    = "Median PFS by CNA Status and Treatment Class",
    subtitle = "Unadjusted — interpret with caution (see adjusted Cox models)\nGroups with n<3 excluded",
    x        = "Treatment Class",
    y        = "Median PFS (months)",
    fill     = "CNA Status"
  ) +
  theme_minimal(base_size=12) +
  theme(
    axis.text.x     = element_text(angle=40, hjust=1),
    strip.text      = element_text(face="bold"),
    legend.position = "top"
  )

ggsave("Figure_C_Median_PFS_by_treatment.png",
       width=10, height=9, dpi=300)
cat("Saved: Figure_C_Median_PFS_by_treatment.png\n")

# ── 14. FIGURE D: KM — CCND1 IN CDK4/6I ONLY ────────────────
cdk_only <- master %>%
  filter(Treatment_Group == "CDK4/6", !is.na(pfs_months))

km_fit <- survfit(
  Surv(pfs_months, pfs_event) ~ CCND1_cna,
  data = cdk_only
)

km_plot <- ggsurvplot(
  km_fit,
  data         = cdk_only,
  pval         = TRUE,
  risk.table   = TRUE,
  title        = "PFS by CCND1 Amplification — CDK4/6 Inhibitor Patients Only",
  legend.labs  = c("CCND1 Wildtype", "CCND1 Amplified"),
  palette      = c("#7FB37F", "#E8836B"),
  xlab         = "Time (months)",
  ylab         = "Progression-free survival probability",
  risk.table.height = 0.25,
  ggtheme      = theme_minimal()
)

ggsave("Figure_D_KM_CCND1_CDK46only.png",
       plot   = print(km_plot),
       width  = 8,
       height = 6,
       dpi    = 300)
cat("Saved: Figure_D_KM_CCND1_CDK46only.png\n")

# ── 15. SAVE ALL RESULTS ─────────────────────────────────────
write_xlsx(
  list(
    "1_Cohort_Overview"        = cohort_overview,
    "2_Univariate_Cox"         = uni_results,
    "3_CCND1_Model"            = ccnd1_results,
    "4_CDKN2A_Model"           = cdkn2a_results,
    "5_EGFR_Model"             = egfr_results,
    "6_Combined_HR"            = combined_hr,
    "7_BH_Correction"          = all_p,
    "8_Descriptive_Summary"    = cna_summary
  ),
  path = "CNA_Analysis_Results_FINAL.xlsx"
)

saveRDS(master, "master_CNA_only_FINAL.rds")
save.image("CNA_analysis_FINAL.RData")

cat("\n============================================================\n")
cat("ALL OUTPUTS SAVED\n")
cat("Excel:   CNA_Analysis_Results_FINAL.xlsx\n")
cat("Figures: Figure_A through Figure_D\n")
cat("RDS:     master_CNA_only_FINAL.rds\n")
cat("============================================================\n")