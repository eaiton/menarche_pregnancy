###############################################################################
#                                                                             #
#            Run main Multivariable Mendelian randomization analyses          #
#                                                                             #
###############################################################################

# Run MVMR of age at menarche on APPOs, with adjustment for adiposity:
# maternal GWAS meta-analysis - Main
# maternal GWAS conditioned on fetal genotype - Sensitivity
# maternal GWAS for single studies - Sensitivity

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
mr_outcome_data_dir <- paste0(home_dir, "/data/b_mr/outcome/")
mvmr_outcome_data_dir <- paste0(home_dir, "/data/c_mvmr/outcome/")
setwd(home_dir)

# Load functions
source("scripts/functions/generate_outdat_with_proxies.R")
source("scripts/functions/run_mvmr.R")

# Select outcomes for each model
outcome_names_df <- read.csv(file.path(mr_outcome_data_dir, "/outcomes.csv"))
# model 1 - all
outcomes_list <- c(outcome_names_df$outcome)
outcome_names_df %>% 
    write.csv(paste0(mvmr_outcome_data_dir, "/outcomes_model_1.csv"),
    row.names = FALSE)

# List of single studies contributing to meta-analysis
studies <- c("ALSPAC", "BIB-SA", "BIB-WE", "FinnGen", "MOBA", "UKB", "Public")

###############################################################################
#                                  Run MVMR                                   #
#               Age at menarche, accounting for body size aged 10             #
###############################################################################

### parameters
# phenotypic covariance matrix to estimate SNP-covariance matrix for
# MVMR package estimates of conditional F statistics
# age at menarche and BMI at age 9 correlation in ALSPAC G1
mvmrcovmatrix <- matrix(c(1, - 0.3012, - 0.3012, 1), nrow = 2, ncol = 2)

### maternal GWAS meta-analysis
outcome_type <- "ma"
outcome_dat_path <- paste0(mr_outcome_data_dir, "/ma_out_dat.txt")
run_mvmr(outcome_type, outcome_dat_path, outcomes_list, mvmrcovmatrix)

### maternal GWAS conditioned on fetal genotype
outcome_type <- "duos"
outcome_dat_path <- paste0(mr_outcome_data_dir, "/duos_out_dat.txt")
run_mvmr(outcome_type, outcome_dat_path, outcomes_list, mvmrcovmatrix)

### maternal GWAS for single studies
for(study_name in studies){
    print(study_name)
    outcome_type <- study_name
    outcome_dat_path <- paste0(mr_outcome_data_dir, "/stu_out_dat_", study_name, ".txt")
    run_mvmr(outcome_type, outcome_dat_path, outcomes_list, mvmrcovmatrix)
}

q('no')