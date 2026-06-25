# HNSCC CNA x Treatment Interaction Analysis

MSc Cancer and Clinical Oncology Dissertation
Barts Cancer Institute / ICR | June 2026

## Summary
Pools per-patient copy number alteration (CNA) and outcome data from
8 HNSCC clinical trials (n=358 patients) to test whether CCND1 amplification,
CDKN2A copy number loss, and EGFR copy number alteration show
treatment-class-specific associations with progression-free survival (PFS).

## Key Findings
- CCND1 amplification: HR=2.89 on CDK4/6 inhibitors (p=0.014, BH p=0.049)
- CDKN2A copy number loss: HR=0.234 on CDK4/6 inhibitors (p=0.002, BH p=0.011)
- CDKN2A x Immunotherapy: HR=6.28 — direction reversal (p=0.002, BH p=0.011)
- EGFR CNA: null result across all treatment classes

## Methods
- CNA encoding: Amp/Del only (mutations excluded)
- Reference group: CDK4/6 inhibitor
- Covariates: HPV status (harmonised)
- Multiple testing: Benjamini-Hochberg FDR correction
- R 4.6.0

## Files
- CNAfinal_analysis.R — complete analysis script
- CNA_Clean_Master (2) copy.xlsx — master dataset (358 patients, 8 studies)
- CNA_Analysis_Results_FINAL.xlsx — all results tables
- Figure_A to Figure_D — key figures
