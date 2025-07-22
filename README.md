# Age at Menarche and Adverse Pregnancy and Perinatal Outcomes: triangulating evidence from multivariable and Mendelian randomization analyses

## Abstract
**Objective** To estimate potential causal effects of age at menarche (AAM) on adverse pregnancy and perinatal outcomes (APPOs).

**Design** Prospective cohort and two-sample Mendelian randomization (MR).

**Setting** European cohorts.

**Population** Up to 934,566 pregnancies from cohorts of predominantly European ancestry with summary genome-wide association data (GWAS) in the MR-PREG collaboration.

**Methods** We estimated confounder adjusted associations using multivariable regression in the Avon Longitudinal Study of Parents and Children (N = 9,441). We used inverse variance weighted analyses for the univariable MR (MR IVW). Sensitivity analyses were performed to test violations of assumptions, including multivariable MR (MVMR) accounting for adiposity.

**Main Outcome Measures** Thirteen APPOs.

**Results** Older AAM was associated with lower risks of hypertensive disorders of pregnancy (HDP), gestational hypertension, preeclampsia, and gestational diabetes mellitus, but adjusting for adiposity attenuated these effects. For example, per 1-year older AAM in MR IVW the OR for HDP was 0.88 (95%CI:0.84, 0.93) and in MVMR was 0.95 (0.90, 1.01), while in multivariable regression the association attenuated from OR=0.91 per 1-year increase in AAM (0.87, 0.94) to OR=0.97 (0.93, 1.01).
We found no clear evidence for effects of AAM on SGA, low birthweight, post-term birth, or perinatal depression from either approach. For other outcomes evidence was limited due to imprecise estimates (very preterm birth), or inconsistent effects in sensitivity analyses (birthweight, large-for-gestational-age, high birthweight, preterm birth).

**Conclusions** We find little robust evidence for causal effects of AAM on APPOs. Effects of younger AAM on increased risks of HDP and GDM may be driven by adiposity.

## Outline
This repo is organised into the following folders and files:
* a_regression: to run multivariable regression analyses
* b_mr: to run univariable Mendelian randomization analyses (MR IVW)
* c_mvmr: to run multivariable Mendelian randomization analyses (MVMR)
* d_combined: to run sensitivities and analyses requiring information across several analyses, including generating figures summarising both approaches
* functions: contains functions required for Mendelian randomization analyses
* instruments: contains lists of SNPs used as genetic instruments in Mendelian randomization analyses
* protocol_24_04_24.pdf: pre-specified analysis plan

## Environmental variables
The code uses some environmental variables which need to be set in your linux environment.
These can be temporarily set with:

```bash
export AGE_AT_MENARCHE_DIR=".../" # directory with /working/**this repo** inside

export MRPREG_sumdat=".../" # directory with outcome GWAS summary statistics saved for extraction
```
