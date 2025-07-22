###############################################################################
#                                                                             #
#                             Run MR IVW analyses                             #
#                                                                             #
###############################################################################

###############################################################################
#                                 Set up                                      #
###############################################################################

# Clear the work environment
rm(list = ls())

# Setting digits
options(digits = 10)

# Required libraries
# remotes::install_github('MRCIEU/TwoSampleMR')
# remotes::install_github("mrcieu/ieugwasr")
x <- c("dplyr", "purrr", "data.table", "TwoSampleMR", "ieugwasr", "metafor",
  "tibble", "broom", "tidyr", "ggplot2")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "/working/")
data_dir <- paste0(home_dir, "data/b_mr/")
gwas_dir <- file.path(data_dir, "gwas/")
exp_data_dir <- file.path(data_dir, "exposure_dat/")
outcome_data_dir <- file.path(data_dir, "outcome/")
formatted_data_dir <- file.path(data_dir, "formatted/")
out_dir <- file.path(home_dir, "results/b_mr/")
scatter_dir <- file.path(out_dir, "plots")
setwd(home_dir)

# Load functions
source("scripts/functions/generate_outdat_with_proxies.R")

###############################################################################
#                        Load exposures & outcomes                            #
###############################################################################

# List of exposures
opengwas_list <- read.csv(paste0(gwas_dir, "/opengwas_list.csv"))
all_opengwasid <- c(opengwas_list$opengwas_id)
all_exposures <- c(opengwas_list$phenotype_id)

# Record pre-selected outcomes
# Abbreviated names used in MR-PREG scripts used to generate variables
outcome_names_df <- read.csv(file.path(outcome_data_dir, "/outcomes.csv"))
all_outcomes <- outcome_names_df$outcome

###############################################################################
#                           Read in outcome datasets                          #
###############################################################################

all_snps <- readLines(file.path(data_dir, "exposure_dat/rsid_all.txt"))

# Maternal GWAS meta-analysis
ma_outdat_all <- read_outcome_data(
            filename = paste0(outcome_data_dir, "/ma_out_dat.txt"),
            snps = all_snps)
# rename sga & lga
ma_outdat_all$outcome[ma_outdat_all$outcome == "lga"] <- "lga_all"
ma_outdat_all$outcome[ma_outdat_all$outcome == "sga"] <- "sga_all"
head(ma_outdat_all)

# Maternal GWAS adjusted for fetal genotype
duos_outdat_all <- fread(paste0(outcome_data_dir, "/duos_out_dat.txt")) %>%
  mutate(chr.outcome = chr, pos.outcome = pos,
  effect_allele.outcome = effect_allele, other_allele.outcome = other_allele,
  eaf.outcome = NA, # to be added in new meta-analysis
  beta.outcome = beta_mat_donuts,
  se.outcome = se_mat_donuts, pval.outcome = p_mat_donuts,
  samplesize.outcome = n_mat, ncase.outcome = NA,
  ncontrol.outcome = NA, outcome = Phenotype, mr_keep.outcome = TRUE,
  pval_origin.outcome = "reported", id.outcome = Phenotype,
  data_source.outcome = "textfile") %>%
  # select columns in ma_outdat
  select(c(colnames(ma_outdat_all))) %>% 
  # for consistency with ma_outdat object
  as.data.frame()
# rename sga & lga
duos_outdat_all$outcome[duos_outdat_all$outcome == "lga"] <- "lga_all"
duos_outdat_all$outcome[duos_outdat_all$outcome == "sga"] <- "sga_all"
head(duos_outdat_all)

# Study-specific data
stu_outdat_all <- fread(paste0(outcome_data_dir, "/stu_out_dat.txt")) %>%
  mutate(chr.outcome = chr, pos.outcome = pos,
  effect_allele.outcome = effect_allele, other_allele.outcome = other_allele,
  eaf.outcome = eaf, beta.outcome = beta, se.outcome = se,
  pval.outcome = pval, samplesize.outcome = samplesize, ncase.outcome = ncase,
  ncontrol.outcome = ncontrol, outcome = Phenotype, mr_keep.outcome = TRUE,
  pval_origin.outcome = "reported", id.outcome = Phenotype,
  data_source.outcome = "textfile") %>%
  # filter to only SNPs in exposure df
  filter(SNP %in% all_snps) %>%
  # remove units in outcome names for consistency
  mutate(outcome = gsub("\\s*\\([^\\)]+\\)","", outcome)) %>%
  # select columns as ma, plus study name
  select(c(colnames(ma_outdat_all), "study")) %>% 
  # for consistency with ma_outdat object
  as.data.frame()
# rename sga & lga
stu_outdat_all$outcome[stu_outdat_all$outcome == "lga"] <- "lga_all"
stu_outdat_all$outcome[stu_outdat_all$outcome == "sga"] <- "sga_all"
head(stu_outdat_all)

# List of included studies
studies <- unique(stu_outdat_all$study)
print(studies)

###############################################################################
#                          Estimate IV strength                               #
###############################################################################

# Using all clumped instruments
# proxies not included since these vary by outcome & dataset

exp_info <- data.frame()
samplesizes_df = data.frame(
  id.exposure = all_exposures[1:2],
  samplesize.exposure = c(632955, 453169)
)

for(exposure_name in all_exposures[1:2]){

    # Read in instruments
    exp_dat <- read.csv(paste0(exp_data_dir, exposure_name,
        "_pval_5e_08_clumped.csv")) %>%
    # Append sample size
      select(-samplesize.exposure) %>%
      left_join(samplesizes_df) %>%
      mutate(
    # Calculate strength
    # Using continuous exposure calculations from Yarmonlinsky et al. 2018
    # PMC6136927
      k = 1,
      maf = pmin(eaf.exposure, 1 - eaf.exposure),
      r2 = 2*(beta.exposure^2)*maf*(1-maf) / (2*(beta.exposure^2)*maf*(1-maf)+se.exposure^2*2*samplesize.exposure*maf*(1-maf)),
      f = r2*(samplesize.exposure-1-k)/((1-r2)*k)
      )

    tmp_exp_info <- group_by(exp_dat, id.exposure) %>%
      summarise(.,
      id    =  unique(id.exposure),
      Nsnps = n(), 
      R2    =  round(sum(r2), 3),
      F     =  round(mean(f), 0),
      max_F = round(max(f), 0),
      min_F = round(min(f), 0),
      median_F = round(median(f), 0),
      )

    exp_info <- rbind(exp_info, tmp_exp_info)
}

write.csv(exp_info, paste0(out_dir, "exp_instrument_info.csv"),
    row.names = FALSE)

###############################################################################
#                              Run MR Analyses                                #
###############################################################################

### 1 - Meta-analysis, including:
###     IVW - Main analysis
###     MR-Egger, Weighted-median, Weighted-mode - Sensitivities
### 2 - Maternal GWAS adjusted for fetal genotype - Sensitivity
### 3 - Single study GWAS - Sensitivity 
### 4 - Leave one study out (LOO) of meta-analysis - Sensitivity
### 5 - Leave out outcome GWAS which overlap with exposure GWAS - Sensitivity

# Set MR methods
mr_methods <-
  c("mr_ivw", # Main analysis
  "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode") # Sensitivities

# Overlapping studies in exposure & outcome GWAS
overlap_studies <- c("UKB", "ALSPAC")

# List of proxy SNPs identified for all exposure SNPs
proxies <- read.csv(paste0(data_dir, "/exposure_dat/proxies.csv"))

# Variables to write to within loop
mr_res <- data.frame()
q_stats <- data.frame()
egger_intercept <- data.frame()
study_loo_res <- data.frame()

# Read in exposure data
exp_dat <- read.csv(paste0(data_dir, "/exposure_dat/aam_pval_5e_08_clumped.csv"))
exposure_name <- "aam"

# Iterate through each outcome, i.e. each exposure-outcome pair
for(outcome_name in all_outcomes){

  print(outcome_name)
  
  # Focussing on each outcome & dataset separately since
  # available proxy SNPs may vary
          
  ### 1 - Meta-analysis - Primary IVW and sensitivities
  
  # Retrieve outcome data, adding proxies where available
  ma_outdat_tmp <- generate_outdat_with_proxies(exp_dat, ma_outdat_all,
    outcome_name, all_outcomes, proxies)
  # Harmonise
  ma_dat <- harmonise_data(exp_dat, ma_outdat_tmp)
  # Run MR
  mr_res_tmp <- mr(ma_dat, method_list = mr_methods) %>% 
    mutate(analysis = "Main", exposure = exposure_name,
    outcome = outcome_name, study = "ma")
  mr_res <- rbind(mr_res, mr_res_tmp)
  # Scatter plots
  p1 <- mr_scatter_plot(mr_res, ma_dat)
  ggsave(p1[[1]], file = paste0(scatter_dir, "/", outcome_name, ".pdf"),
    width = 7, height = 7)
  # Cochrane's Q statistic
  q_stats_tmp <- mr_heterogeneity(ma_dat) %>% 
    mutate(analysis = "Main", exposure = exposure_name,
    outcome = outcome_name, study = "ma")
  q_stats <- rbind(q_stats, q_stats_tmp)
  # MR-Egger intercept test
  egger_intercept_tmp <- mr_pleiotropy_test(ma_dat) %>% 
    mutate(analysis = "Main", exposure = exposure_name,
    outcome = outcome_name, study = "ma")
  egger_intercept <- rbind(egger_intercept, egger_intercept_tmp)
  
  ### 2 - Maternal GWAS adjusted for fetal genotype - Sensitivity
  
  # Retrieve outcome data, adding proxies where available
  duos_outdat_tmp <- generate_outdat_with_proxies(exp_dat, duos_outdat_all,
    outcome_name, all_outcomes, proxies)
  if(dim(duos_outdat_tmp)[1] == 0){
    print("Outcome not available in fetal genotype adjusted GWAS")
  } else {
    # Harmonise
    duos_dat <- harmonise_data(exp_dat, duos_outdat_tmp)
    # Run MR
    duos_res_tmp <- mr(duos_dat, method_list = mr_methods) %>% 
      mutate(analysis = "Adjusted for fetal genotype",
      exposure = exposure_name, outcome = outcome_name, study = "duos")
    mr_res <- rbind(mr_res, duos_res_tmp)
    # Cochrane's Q statistic
    duos_q_stats_tmp <- mr_heterogeneity(duos_dat) %>% 
      mutate(analysis = "Adjusted for fetal genotype",
      exposure = exposure_name, outcome = outcome_name, study = "duos")
    q_stats <- rbind(q_stats, duos_q_stats_tmp)
    # MR-Egger intercept test
    duos_egger_intercept_tmp <- mr_pleiotropy_test(duos_dat) %>% 
      mutate(analysis = "Adjusted for fetal genotype",
      exposure = exposure_name, outcome = outcome_name, study = "duos")
    egger_intercept <- rbind(egger_intercept, duos_egger_intercept_tmp)
  }
  
  ### 3 - Single study GWAS - Sensitivity
  
  stu_res <- data.frame()
  
  for(study_name in studies){
  
    # Filter outcome data to study, remove 'study' column
    stu_outdat_all_tmp <- stu_outdat_all %>% 
      filter(study == study_name) %>% 
      select(-study)
    # Retrieve outcome data, adding proxies where available
    stu_outdat_tmp <- generate_outdat_with_proxies(exp_dat, stu_outdat_all_tmp,
      outcome_name, all_outcomes, proxies)
    if(dim(stu_outdat_tmp)[1] == 0){
      print(paste("Outcome", outcome_name, "not available in", study_name))
    } else if(dim(stu_outdat_tmp)[1] > 0){
      # Harmonise
      stu_outdat_tmp$id.outcome <- stu_outdat_tmp$id.outcome[1]
      stu_dat <- harmonise_data(exp_dat, stu_outdat_tmp)
      # Run MR
      stu_res_tmp <- mr(stu_dat, method_list = mr_methods) %>% 
        mutate(analysis = paste(study_name, "only"), exposure = exposure_name,
        outcome = outcome_name, study = study_name)
      # Save each study specific result for this exposure-outcome pair,
      # to be used in subsequent sensitivities
      stu_res <- rbind(stu_res, stu_res_tmp)
    }

}

mr_res <- rbind(mr_res, stu_res)

### 4 - Leave one study out (LOO) of meta-analysis - Sensitivity

# Format data for LOO
stu_res_ivw <- filter(stu_res, method == "Inverse variance weighted")
study_loo_dat <- metafor::rma(yi = b, sei = se, data = stu_res_ivw,
 slab = study, method = "FE")
# Run LOO removing one study at a time
print("Running leave one study out (LOO) MR analysis")
study_loo_res_tmp <- metafor::leave1out(study_loo_dat)
# Format LOO output
study_loo_res_tmp_formatted <- study_loo_res_tmp %>% 
  as.data.frame %>%
  tibble::rownames_to_column(., var = "study") %>%
  mutate(se = estimate/zval,
  analysis = paste(study, "left out"),
  exposure = exposure_name,
  outcome = outcome_name)
study_loo_res <- rbind(study_loo_res, study_loo_res_tmp_formatted)

### 5 - Leave out overlap GWAS - Sensitivity
# Format 
no_overlap_tmp <- filter(stu_res_ivw, !(study %in% overlap_studies))
# Meta-analysis
no_overlap_res_tmp <- metafor::rma(yi = b, sei = se,
  data = no_overlap_tmp, method = "FE")
no_overlap_res_tmp_2 <- no_overlap_res_tmp %>%
  tidy() %>% as.data.frame %>%
# to match first loo analysis columns
  mutate(analysis = "UKB & ALSPAC left out",
    exposure = exposure_name, outcome = outcome_name,
    se = std.error, zval = NA, pval = p.value,
    study = "No overlap") 
no_overlap_res_tmp_2$ci.lb <- no_overlap_res_tmp$ci.lb
no_overlap_res_tmp_2$ci.ub <- no_overlap_res_tmp$ci.ub
no_overlap_res_tmp_2$Q <- no_overlap_res_tmp$QE
  # test statistic of the test for (residual) heterogeneity
no_overlap_res_tmp_2$Qp <- no_overlap_res_tmp$QEp
no_overlap_res_tmp_2$I2 <- no_overlap_res_tmp$I2
no_overlap_res_tmp_2$H2 <- no_overlap_res_tmp$H2
no_overlap_res_tmp_2$se <- no_overlap_res_tmp$se
# save same columns as LOO analysis
no_overlap_res_tmp_2 <- no_overlap_res_tmp_2 %>% 
  select(colnames(study_loo_res))
  study_loo_res <- rbind(study_loo_res, no_overlap_res_tmp_2)

}

###############################################################################
#                         Append sample sizes                                 #
###############################################################################

### Record sample sizes for outcome & datasets
### Using medians since sample size varies for each SNP instrument
### Should be re-run once UKB sample size added

### Main
ma_sample_sizes <- ma_outdat_all %>%
  group_by(outcome) %>%
  summarise(min_n = min(samplesize.outcome),
  mean_n = mean(samplesize.outcome),
  median_n = median(samplesize.outcome),
  max_n = max(samplesize.outcome),
  median_ncases = median(ncase.outcome),
  median_ncontrol = median(ncontrol.outcome)) %>% 
  # round to nearest integers
  mutate_if(is.numeric, round) %>% 
  mutate(analysis = "Main")

### Adjusted for fetal genotype
duos_sample_sizes <- duos_outdat_all %>%
  group_by(outcome) %>%
  summarise(min_n = min(samplesize.outcome),
  mean_n = mean(samplesize.outcome),
  median_n = median(samplesize.outcome),
  max_n = max(samplesize.outcome),
  # case & control data not available for this dataset,
  # but save NAs to allow rbind()
  median_ncases = NA,
  median_ncontrol = NA) %>% 
  # round to nearest integers
  mutate_if(is.numeric, round) %>% 
  mutate(analysis = "Adjusted for fetal genotype")

### Study-level individual analysis
stu_sample_sizes <- stu_outdat_all %>%
  group_by(outcome, study) %>%
  summarise(min_n = min(samplesize.outcome),
  mean_n = mean(samplesize.outcome),
  median_n = median(samplesize.outcome),
  max_n = max(samplesize.outcome),
  median_ncases = median(ncase.outcome),
  median_ncontrol = median(ncontrol.outcome)) %>% 
  # round to nearest integers
  mutate_if(is.numeric, round) %>% 
  mutate(analysis = paste(study, "only")) %>%
  select(-study)

mr_sample_sizes <- rbind(ma_sample_sizes, duos_sample_sizes, stu_sample_sizes)

### Leave one study out meta-analysis
calculate_loo_sample_size <- function(study){
  if(study == "UKB & ALSPAC"){
    tmp <- stu_sample_sizes %>% 
    filter(analysis != "UKB only" & analysis != "ALSPAC only")
  } else {
    tmp <- stu_sample_sizes %>% 
    filter(analysis != paste(study, "only"))
  }
  tmp %>%
    group_by(outcome) %>%
    summarise(sum_min_n = sum(min_n, na.rm = TRUE),
    sum_mean_n = sum(mean_n, na.rm = TRUE),
    sum_median_n = sum(median_n, na.rm = TRUE),
    sum_max_n = sum(max_n, na.rm = TRUE),
    sum_median_ncases = sum(median_ncases, na.rm = TRUE),
    sum_median_ncontrol = sum(median_ncontrol, na.rm = TRUE)) %>%
    # round to nearest integers
    mutate_if(is.numeric, round) %>%
    mutate(analysis = paste(study, "left out"))
}

loo_sample_sizes <- data.frame()
for(study_name in c(studies, "UKB & ALSPAC")){ 
  loo_tmp <- calculate_loo_sample_size(study_name)
  loo_sample_sizes <- rbind(loo_sample_sizes, loo_tmp)
}

###############################################################################
#                         Write out all MR results                            #
###############################################################################

### Write out results for all exposure & outcome pairs tested,
### appending full outcome names

mr_res %>%
  left_join(outcome_names_df) %>%
  left_join(mr_sample_sizes, by = join_by(outcome, analysis)) %>%
  write.csv(paste0(out_dir, "/mr_results.csv"), row.names = FALSE)

q_stats %>% 
  left_join(outcome_names_df) %>% 
  write.csv(paste0(out_dir, "/q_statistics.csv"), row.names = FALSE)

egger_intercept %>% 
  left_join(outcome_names_df) %>%
  write.csv(paste0(out_dir, "/egger_intercept.csv"), row.names = FALSE)

study_loo_res %>% 
  left_join(outcome_names_df) %>%
  left_join(loo_sample_sizes, by = join_by(outcome, analysis)) %>%
  write.csv(paste0(out_dir, "/leave_one_study_out_results.csv"),
    row.names = FALSE)

### Add columns for plots & tables - separate binary & continuous
mr_res %>%
  left_join(outcome_names_df) %>%
  left_join(mr_sample_sizes, by = join_by(outcome, analysis)) %>%
  # remove birthweight
  filter(outcome != "zbw_all") %>%
  # add CIs
  generate_odds_ratios() %>%
  mutate(
    approach = "MR Univariate",
    `Cases / controls` = paste0(median_ncases, " / ", median_ncontrol),
    # generate clean OR label
    `OR (95% CI)` = sprintf("%.2f (%.2f - %.2f)", or, or_lci95, or_uci95)
   ) %>%
  write.csv(paste0(out_dir, "/mr_results_binary.csv"), row.names = FALSE)

mr_res %>%
  left_join(outcome_names_df) %>%
  left_join(mr_sample_sizes, by = join_by(outcome, analysis)) %>%
  # restrict to birthweight
  filter(outcome == "zbw_all") %>%
  # add CIs
  mutate(
    approach = "MR Univariate",
    lo_ci = b - qnorm(0.975)*se,
    up_ci = b + qnorm(0.975)*se,
    `Sample size` = paste0(median_n),
    # generate clean beta label
     `Beta (95% CI)` = sprintf("%.2f (%.2f - %.2f)", b, lo_ci, up_ci)
     ) %>%
  write.csv(paste0(out_dir, "/mr_results_bw.csv"), row.names = FALSE)

q('no')