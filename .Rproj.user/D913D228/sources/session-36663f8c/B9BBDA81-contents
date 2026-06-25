# ============================================================
# CNA × Treatment Interaction Analysis — HNSCC
# MSc Cancer and Clinical Oncology Dissertation
# ============================================================

# ── 0. LIBRARIES ─────────────────────────────────────────────
library(readxl)
library(dplyr)
library(tidyr)
library(survival)
library(broom)
library(ggplot2)
library(writexl)
library(survminer)
library(tibble)

# ── 1. LOAD DATA ─────────────────────────────────────────────
master_raw <- read_excel(
  "CNA_Clean_Master (2) copy.xlsx",
  sheet = "Master Matrix",
  col_names = FALSE
)

colnames(master_raw) <- as.character(master_raw[2, ])
master <- master_raw[-c(1, 2), ]

# ── 2. CLEAN COLUMN NAMES ────────────────────────────────────
master <- master %>%
  rename(
    patient_id = `Patient ID`,
    study = `Study`,
    cancer_type = `Cancer\r\nType`,
    treatment_class = `Treatment Class`,
    drug = `Drug`,
    setting = `Setting`,
    hpv = `HPV`,
    bor = `BOR`,
    clinical_benefit = `Clinical\r\nBenefit`,
    pfs_months = `PFS\r\n(mo)`,
    os_months = `OS\r\n(mo)`
  ) %>%
  mutate(
    pfs_months = as.numeric(pfs_months),
    os_months = as.numeric(os_months)
  )

# ── 3. CNA BINARY VARIABLES ──────────────────────────────────
master <- master %>%
  mutate(
    CCND1_cna  = ifelse(CCND1 == "Amp", 1, 0),
    CDKN2A_cna = ifelse(CDKN2A %in% c("Amp", "Del"), 1, 0),
    EGFR_cna   = ifelse(EGFR %in% c("Amp", "Del"), 1, 0),
    PIK3CA_cna = ifelse(PIK3CA == "Amp", 1, 0)
  )

# ── 4. TREATMENT GROUPS ──────────────────────────────────────
master <- master %>%
  mutate(
    Treatment_Group = case_when(
      grepl("CDK4/6", treatment_class) ~ "CDK4/6",
      grepl("EGFR", treatment_class) ~ "EGFR",
      grepl("Immunotherapy", treatment_class) ~ "Immunotherapy",
      grepl("PI3K", treatment_class) ~ "PI3K",
      TRUE ~ "Other"
    )
  )

# HPV clean variable (FIX)
master$HPV_clean <- master$hpv

# ── 5. COHORT OVERVIEW ───────────────────────────────────────
cohort_overview <- master %>%
  group_by(Treatment_Group) %>%
  summarise(
    N = n(),
    PFS_available = sum(!is.na(pfs_months)),
    OS_available = sum(!is.na(os_months)),
    CCND1_CNA_n = sum(CCND1_cna, na.rm=TRUE),
    CDKN2A_CNA_n = sum(CDKN2A_cna, na.rm=TRUE),
    EGFR_CNA_n = sum(EGFR_cna, na.rm=TRUE),
    PIK3CA_CNA_n = sum(PIK3CA_cna, na.rm=TRUE),
    .groups = "drop"
  )

print(cohort_overview)

# ── 6. UNIVARIATE COX MODELS ─────────────────────────────────
genes_cna <- c("CCND1_cna", "CDKN2A_cna", "EGFR_cna", "PIK3CA_cna")

uni_results <- lapply(genes_cna, function(g) {
  fit <- coxph(as.formula(paste0("Surv(pfs_months) ~ ", g)), data = master)
  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(gene = g)
}) %>%
  bind_rows() %>%
  select(gene, estimate, conf.low, conf.high, p.value) %>%
  rename(HR = estimate, CI_low = conf.low, CI_high = conf.high, p = p.value)

print(uni_results)

# ── 7. INTERACTION MODELS ─────────────────────────────────────
master_model <- master %>%
  filter(Treatment_Group != "Other", !is.na(pfs_months))

master_model$Treatment_Group <- relevel(
  factor(master_model$Treatment_Group),
  ref = "CDK4/6"
)

cox_ccnd1 <- coxph(Surv(pfs_months) ~ CCND1_cna * Treatment_Group + HPV_clean,
                   data = master_model)

cox_cdkn2a <- coxph(Surv(pfs_months) ~ CDKN2A_cna * Treatment_Group + HPV_clean,
                    data = master_model)

cox_egfr <- coxph(Surv(pfs_months) ~ EGFR_cna * Treatment_Group + HPV_clean,
                  data = master_model)

ccnd1_results <- tidy(cox_ccnd1, exponentiate=TRUE, conf.int=TRUE)
cdkn2a_results <- tidy(cox_cdkn2a, exponentiate=TRUE, conf.int=TRUE)
egfr_results <- tidy(cox_egfr, exponentiate=TRUE, conf.int=TRUE)

# ── 8. COMBINED HR TABLE ─────────────────────────────────────
combined_hr <- tibble(
  Gene = c("CCND1", "CCND1", "CCND1", "CDKN2A", "CDKN2A"),
  Treatment = c("CDK4/6", "EGFR", "Immunotherapy", "CDK4/6", "Immunotherapy"),
  HR = c(
    ccnd1_results$estimate[1],
    ccnd1_results$estimate[1] * ccnd1_results$estimate[2],
    ccnd1_results$estimate[1] * ccnd1_results$estimate[3],
    cdkn2a_results$estimate[1],
    cdkn2a_results$estimate[1] * cdkn2a_results$estimate[3]
  )
)

print(combined_hr)

# ── 9. BH CORRECTION ─────────────────────────────────────────
all_p <- tibble(
  test = c(
    "CCND1 main", "CDKN2A main", "EGFR main", "PIK3CA main",
    "CCND1 × EGFR", "CCND1 × Immunotherapy",
    "CDKN2A × Immunotherapy", "EGFR × EGFR"
  ),
  p_raw = c(
    uni_results$p[1],
    uni_results$p[2],
    uni_results$p[3],
    uni_results$p[4],
    ccnd1_results$p.value[2],
    ccnd1_results$p.value[3],
    cdkn2a_results$p.value[3],
    egfr_results$p.value[2]
  )
) %>%
  filter(!is.na(p_raw)) %>%
  mutate(
    p_BH = p.adjust(p_raw, method = "BH"),
    sig = p_BH < 0.05
  )

print(all_p)

# ── 10. DESCRIPTIVE CNA TABLE ────────────────────────────────
cna_summary <- master %>%
  filter(Treatment_Group != "Other") %>%
  pivot_longer(cols = c(CCND1_cna, CDKN2A_cna, EGFR_cna),
               names_to = "Gene",
               values_to = "CNA_status") %>%
  mutate(
    Gene = recode(Gene,
                  "CCND1_cna" = "CCND1",
                  "CDKN2A_cna" = "CDKN2A",
                  "EGFR_cna" = "EGFR"),
    CNA_label = ifelse(CNA_status == 1, "CNA+", "CNA−")
  ) %>%
  group_by(Gene, Treatment_Group, CNA_label) %>%
  summarise(
    n = n(),
    CB_n = sum(clinical_benefit == "Benefit", na.rm=TRUE),
    CB_pct = round(100 * CB_n / n, 1),
    median_pfs = median(pfs_months, na.rm=TRUE),
    .groups = "drop"
  )

print(cna_summary)

