###############################################################################
#                     Function to run MVMR using local                        #
#                 proxies and MR-PREG outcome data formats                    #
###############################################################################

# MVMR with 2 exposures: age at menarche & adiposity

# Selects outcome SNPs based on availability in outcome data,
# using proxies when needed

# Generates conditional F stats, Q stats, scatter plots,
# MR effect estimate results from both TwoSampleMR package and MVMR package 

# Requires: library(dplyr), library(tidyr), library(TwoSampleMR), library(MVMR)

run_mvmr <- function(outcome_type, outcome_dat_path, outcomes_list,
  mvmrcovmatrix){

### Make directories to save results ##########################################
out_dir <- paste0(home_dir, "/results/c_mvmr/", outcome_type)
plots_out_dir <- paste0(out_dir, "/plots/")
raw_out_dir <- paste0(out_dir, "/raw/")
system(paste0("mkdir -p ", out_dir))
system(paste0("mkdir -p ", plots_out_dir))
system(paste0("mkdir -p ", raw_out_dir))

### Read in data ##############################################################
### Exposure
exposure_dat <-
  read.csv(paste0(home_dir, "/data/c_mvmr/instruments/mvmr_snps.csv"))
  
### Proxies
# Use proxies file generated by scripts/b_mr/2_b_sel_snps.R
proxies <- read.csv(paste0(home_dir, "/data/b_mr/exposure_dat/proxies.csv"))
proxy_snps <- unique(proxies$rsid)
length(proxy_snps)

### Outcomes
outcome_names_df <- read.csv(
  paste0(home_dir, "/data/c_mvmr/outcome/outcomes.csv"))

### Outcome data
all_snps <- c(exposure_dat$SNP, proxy_snps)

if(grepl("duos", outcome_dat_path) == TRUE){ # fetal-adjusted maternal GWAS
    outcome_dat <- read_outcome_data(
    filename = outcome_dat_path,
    snps = all_snps,
    sep = " ",
    phenotype_col = "Phenotype",
    snp_col = "SNP",
    beta_col = "beta_mat_donuts",
    se_col = "se_mat_donuts",
    eaf_col = "eaf_mat",
    pval_col = "p_mat_donuts",
    samplesize_col = "n_mat"
  )
} else if(grepl("stu", outcome_dat_path) == TRUE){ # single-study maternal GWAS
  outcome_dat <- read_outcome_data(
    filename = outcome_dat_path,
    snps = all_snps,
    sep = " ") %>%
    # remove units in outcome names
    mutate(outcome = stringr::str_replace(outcome, "(?s) .*", ""))
} else { # meta-analysis maternal GWAS
  outcome_dat <- read_outcome_data(
    filename = outcome_dat_path,
    snps = all_snps,
    sep = " ")
}
# rename sga & lga
outcome_dat$outcome[outcome_dat$outcome == "lga"] <- "lga_all"
outcome_dat$outcome[outcome_dat$outcome == "sga"] <- "sga_all"
# keep only outcomes of interest
outcome_dat <- outcome_dat[outcome_dat$outcome %in% outcomes_list, ]
# stop running if no SNPs for outcomes of interest available
if(dim(outcome_dat)[1] == 0){
  warning("No outcomes of interest available for MVMR")
}
# arrange in same order for consistency with observational results
outcome_dat <- outcome_dat %>% arrange(match(outcome, outcomes_list))
# snps avilable before proxies
outcome_dat %>%
  filter(SNP %in% exposure_dat$SNP) %>%
  count(outcome)

### Run MVMR using proxies ####################################################

# Available outcomes of interest
avail_outcomes <- unique(outcome_dat$outcome)

# Empty dfs/lists to save results in
mvmr_results_twosample <- data.frame()
aam_plot_list <- list()
bodysize_plot_list <- list()
f_stats <- data.frame()
q_stats <- data.frame()
mvmr_results_mvmr <- data.frame()
sample_sizes <- data.frame()

for (i in seq_along(avail_outcomes)) {

    # For testing
    i_outcome <- avail_outcomes[i]
    print(i_outcome)

    ### 0. add proxies ########################################################
    tmp_outcome_dat <- generate_outdat_with_proxies(exposure_dat, outcome_dat,
    i_outcome, avail_outcomes, proxies)
    print("proxies added")

	  ### 1. harmonise ##########################################################
	  tmp_outcome_dat$id.outcome <- tmp_outcome_dat$id.outcome[1]
    tmp_outcome_dat$data_source.outcome <- "textfile"
    tmp_mvdat <- mv_harmonise_data(exposure_dat, tmp_outcome_dat)
    print("data harmonised")
	
	  ### 2.a run mvmr with TwoSampleMR package #################################
	  res_two_sample <- mv_multiple(tmp_mvdat, plots = TRUE)
    print("twosamplemr run")
	  # save causal effect estimates
	  mvmr_results_twosample <- rbind(mvmr_results_twosample,
	    res_two_sample$result)
	  # save plots to list
    bodysize_plot_list[[i]] <- res_two_sample$plots[[1]]
	  aam_plot_list[[i]] <- res_two_sample$plots[[2]]
	
	  ### 2.b run mvmr with mvmr package ########################################
	  # re-format summary data using harmonisation from TwoSampleMR::mv_multiple
    # x1 = adiposity
    # x2 = age at menarche
	  F.data <- format_mvmr(
	    BXGs = tmp_mvdat$exposure_beta[, c(1, 2)],
	    BYG = tmp_mvdat$outcome_beta,
	    seBXGs = tmp_mvdat$exposure_se[, c(1, 2)],
	    seBYG = tmp_mvdat$outcome_se,
	    RSID = rownames(tmp_mvdat$exposure_beta))
    print("mvmr data formatted")
	
	  # test for weak instruments - F statistics (Qx)
	  # assuming no covariance
	  sres_no_cov <- strength_mvmr(r_input = F.data, gencov = 0) %>%
	    mutate(adj_for_cov = FALSE,
	            outcome = i_outcome)
	  # using covariance matrix
	  Xcovmat <- phenocov_mvmr(mvmrcovmatrix, seBXGs = F.data[, 6:7])
	  sres_adj_cov <- strength_mvmr(r_input = F.data, gencov = Xcovmat) %>%
	    mutate(adj_for_cov = TRUE,
	           outcome = i_outcome)
	  f_stats <- rbind(f_stats, sres_no_cov, sres_adj_cov)
    print("f stats calculated")
	
	  # test for horizontal pleiotropy - Q statistic (Qa)
    # degrees of freedom printed in .Rout but not saved in object
	  pres <- pleiotropy_mvmr(r_input = F.data, gencov = Xcovmat) %>%
	    as.data.frame() %>%
	   mutate(adj_for_cov = TRUE,
	          outcome = i_outcome)
	  q_stats <- rbind(q_stats, pres)
    print("q stats calculated")
	
	  # run mvmr and save results
	  # assuming no covariance
	  res_mvmr <- ivw_mvmr(r_input = F.data) %>%
	    as.data.frame() %>%
	    mutate(adj_for_cov = FALSE,
	            outcome = i_outcome) %>%
	    tibble::rownames_to_column("exposure")
	  # using covariance matrix
	  res_mvmr_cov <- ivw_mvmr(r_input = F.data, gencov = Xcovmat) %>%
	    as.data.frame() %>%
	    mutate(adj_for_cov = TRUE,
	            outcome = i_outcome) %>%
	    tibble::rownames_to_column("exposure")
	  # save both
	  mvmr_results_mvmr <- rbind(mvmr_results_mvmr,
	    res_mvmr, res_mvmr_cov)
    print("MVMR package run")

    ### 3. save tmp_outcome_dat sample sizes
    ### all data available for single studies & meta-analysis
    ### (except UKB in release 2)
    if(grepl("ma", outcome_dat_path) == TRUE | grepl("stu", outcome_dat_path) == TRUE){
    tmp_sample_sizes <- tmp_outcome_dat %>%
      group_by(outcome) %>%
      summarise(min_n = min(samplesize.outcome, na.rm = TRUE),
                mean_n = mean(samplesize.outcome, na.rm = TRUE),
                median_n = median(samplesize.outcome, na.rm = TRUE),
                max_n = max(samplesize.outcome, na.rm = TRUE),
                median_ncases = median(ncase.outcome, na.rm = TRUE),
                median_ncontrol = median(ncontrol.outcome, na.rm = TRUE)) %>%
      # round to nearest integer
      mutate_if(is.numeric, round)
    sample_sizes <- rbind(sample_sizes, tmp_sample_sizes)
    ### for conditional estimates, no case/control numbers in release 2
    } else if(grepl("duos", outcome_dat_path) == TRUE){
    tmp_sample_sizes <- tmp_outcome_dat %>%
      group_by(outcome) %>%
      summarise(min_n = min(samplesize.outcome, na.rm = TRUE),
                mean_n = mean(samplesize.outcome, na.rm = TRUE),
                median_n = median(samplesize.outcome, na.rm = TRUE),
                max_n = max(samplesize.outcome, na.rm = TRUE),
                median_ncases = NA,
                median_ncontrol = NA) %>%
      # round to nearest integer
      mutate_if(is.numeric, round)
    sample_sizes <- rbind(sample_sizes, tmp_sample_sizes)
    }
    print("sample sizes saved")

}

### Save results ##############################################################

# Save adiposity plots, separate file for each plot
for (i in seq_along(avail_outcomes)) {
  i_outcome <- avail_outcomes[i]
  file_name <- paste0(plots_out_dir, i_outcome, "_body_size.png")
  png(file_name, width = 800, height = 1200)
  print(bodysize_plot_list[[i]])
  dev.off()
}
# Save aam plots, separate file for each plot
for (i in seq_along(avail_outcomes)) {
  i_outcome <- avail_outcomes[i]
  file_name <- paste0(plots_out_dir, i_outcome, "_age_at_menarche.png")
  png(file_name, width = 800, height = 1200)
  print(aam_plot_list[[i]])
  dev.off()
}
print("saved plots")

if(dim(outcome_dat)[1] > 0) {

# Save numeric results
head(mvmr_results_twosample)
mvmr_results_twosample %>%
  left_join(outcome_names_df, by = join_by(outcome)) %>%
  left_join(sample_sizes, by = join_by(outcome)) %>%
  write.csv(paste0(raw_out_dir, "/mvmr_results_twosample.csv"),
    row.names = FALSE)

f_stats %>%
  left_join(outcome_names_df, by = join_by(outcome)) %>%
  write.csv(paste0(raw_out_dir, "/f_stats_mvmr.csv"),
    row.names = FALSE)

q_stats %>% 
  left_join(outcome_names_df, by = join_by(outcome)) %>%
  mutate(
    Qstat = signif(Qstat, 3),
    Qpval = signif(Qpval, 3)) %>%
  select(outcome_clean, Qstat, Qpval) %>% 
  write.csv(paste0(raw_out_dir, "/q_stats_mvmr.csv"),
          row.names = FALSE)

### Change exposure names on MVMR package (match covariance matrix)
head(mvmr_results_mvmr)
mvmr_results_mvmr <- mvmr_results_mvmr %>%
  mutate(
    exposure = case_when(exposure == "exposure1" ~ "body_size_10",
    exposure == "exposure2" ~ "aam"
    )
  )
head(mvmr_results_mvmr)
write.csv(mvmr_results_mvmr,
          paste0(raw_out_dir, "/mvmr_results_mvmr.csv"),
          row.names = FALSE)

### checking consistency results across both packages:
mvmr_results_twosample[mvmr_results_twosample$outcome == "zbw_all", ]
mvmr_results_mvmr[mvmr_results_mvmr$outcome == "zbw_all", ]

} else { warning("No MVMR results saved") }

print("saved results")

}