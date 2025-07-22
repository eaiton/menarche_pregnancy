###############################################################################
#                                                                             #
#                                Clean results from                           #
#                     Multivariable Mendelian randomization                   #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list=ls())

# Setting digits
options(digits = 10)

# Required libraries
## install.packages("remotes")
## remotes::install_github("MRCIEU/TwoSampleMR")
## remotes::install_github("WSpiller/MVMR")
x <- c("dplyr", "tidyr", "TwoSampleMR", "MVMR", "stringr", "metafor", "broom")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "/working/")
#mr_results_dir <- paste0(home_dir, "/results/mr/outcome/")
mvmr_results_dir <- paste0(home_dir, "/results/c_mvmr/")
setwd(home_dir)

###############################################################################
#                           Combine all results                               #
###############################################################################

mvmr_results_folders <- list.dirs(mvmr_results_dir)

### For all main results, like MR
mvmr_results_folders_main <-
    mvmr_results_folders[grepl("raw$", mvmr_results_folders) &
    !grepl("archive", mvmr_results_folders) & !grepl("snp", mvmr_results_folders)]
    # remove snp-test from GDM issue testing
mvmr_results_folders_main

mvmr_res <- data.frame()
f_stats <- data.frame()
q_stats <- data.frame()

for(folder in mvmr_results_folders_main){

    print(folder)

    if(grepl("ma", folder)){
        analysis_name = "Main"
    } else if (grepl("duos", folder)){
        analysis_name = "Adjusted for fetal genotype"
    } else {
        analysis_name = paste(sub(".*/(.*)/raw$", "\\1", folder), "only")
    }
    
    # read in results and append info
    res_tmp <- read.csv(paste0(folder, "/mvmr_results_twosample.csv")) %>%
        mutate(analysis = analysis_name)
    f_stats_tmp <- read.csv(paste0(folder, "/f_stats_mvmr.csv")) %>%
        mutate(analysis = analysis_name)
    if(model_name == "model_1") { f_stats_tmp$exposure3 <- NA }
    q_stats_tmp <- read.csv(paste0(folder, "/q_stats_mvmr.csv")) %>%
        mutate(analysis = analysis_name)

    mvmr_res <- rbind(mvmr_res, res_tmp)
    f_stats <- rbind(f_stats, f_stats_tmp)
    q_stats <- rbind(q_stats, q_stats_tmp)
    
}

### add leave one study out results
mvmr_res_loo <- read.csv(
    paste0(mvmr_results_dir, "/model_1/loo/leave_one_study_out_results.csv")) %>%
    mutate(model = "model_1", b = estimate, up_ci = ci.ub, lo_ci = ci.lb)

# to allow for rbind() - add all columns to both
cols_to_add_to_main <- colnames(mvmr_res_loo)[!(colnames(mvmr_res_loo) %in% colnames(mvmr_res))]
for(col in cols_to_add_to_main){
    mvmr_res[[col]] <- NA
}
cols_to_add_to_loo <- colnames(mvmr_res)[!(colnames(mvmr_res) %in% colnames(mvmr_res_loo))]
for(col in cols_to_add_to_loo){
    mvmr_res_loo[[col]] <- NA
}

mvmr_res <- rbind(mvmr_res, mvmr_res_loo)

### Save these combined files
mvmr_res %>%
    write.csv(paste0(mvmr_results_dir, "/mvmr_results.csv"),
    row.names = FALSE)
f_stats %>%
    write.csv(paste0(mvmr_results_dir, "/f_statistics.csv"),
    row.names = FALSE)
q_stats %>%
    write.csv(paste0(mvmr_results_dir, "/q_statistics.csv"),
    row.names = FALSE)

### Add columns for plots & tables - separate binary & continuous
mvmr_res %>%
  # remove birthweight
  filter(outcome != "zbw_all") %>%
  # add CIs
  generate_odds_ratios() %>%
  mutate(
    approach = "MR Multivariable",
    `Cases / controls` = paste0(median_ncases, " / ", median_ncontrol),
    # generate clean OR label
    `OR (95% CI)` = sprintf("%.2f (%.2f - %.2f)", or, or_lci95, or_uci95)
   ) %>%
  write.csv(paste0(mvmr_results_dir, "/mvmr_results_binary.csv"), row.names = FALSE)

mvmr_res %>%
  # restrict to birthweight
  filter(outcome == "zbw_all") %>%
  # add CIs
  mutate(
    approach = "MR Multivariable",
    lo_ci = b - qnorm(0.975)*se,
    up_ci = b + qnorm(0.975)*se,
    `Sample size` = paste0(median_n),
    # generate clean beta label
     `Beta (95% CI)` = sprintf("%.2f (%.2f - %.2f)", b, lo_ci, up_ci)
     ) %>%
  write.csv(paste0(mvmr_results_dir, "/mvmr_results_bw.csv"), row.names = FALSE)

q('no')