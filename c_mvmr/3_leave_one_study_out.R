###############################################################################
#                                                                             #
#                     Meta-analysis sensitivities for                         #
#              Multivarible Mendelian randomization analyses                  #
#                                                                             #
###############################################################################

# Run MVMR of age at menarche on APPOs, with adjustment for adiposity:
# leave one study out - Sensitivity
# leave out overlapping studies in both exposure & outcome GWAS  - Sensitivity

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list=ls())

# Setting digits
options(digits = 10)

# Required libraries
x <- c("dplyr", "tidyr", "stringr", "metafor", "broom")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "/working/")
mr_outcome_data_dir <- paste0(home_dir, "/data/b_mr/outcome/")
mvmr_outcome_data_dir <- paste0(home_dir, "/data/c_mvmr/outcome/")
setwd(home_dir)

# Select outcomes for each model
outcome_names_df <- read.csv(file.path(mr_outcome_data_dir, "/outcomes.csv"))
# model 1 - all
model_1_outcomes_list <- c(outcome_names_df$outcome)

# List of single studies contributing to meta-analysis
ma_out <- read.table(paste0(mr_outcome_data_dir, "/stu_out_dat.txt"), header = TRUE)
studies <- unique(ma_out$study)
rm(ma_out) # to save memory

###############################################################################
#                       Sample size calculation set up                        #
###############################################################################

# Overlapping studies in exposure & outcome GWAS
overlap_studies <- c("UKB", "ALSPAC")

# Function to calculate sample sizes for MVMR leave one study out meta-analysis
# study_name = the study being left out
calculate_loo_sample_size <- function(study_name){
  if(study_name == "UKB & ALSPAC"){
    tmp <- stu_res_all %>% 
    filter(!(analysis %in% c("UKB", "ALSPAC")))
  } else {
    tmp <- stu_res_all %>% 
    filter(analysis != study_name)
  }
  tmp %>%
    group_by(outcome, id.exposure) %>%
    summarise(sum_min_n = sum(min_n, na.rm = TRUE),
    sum_mean_n = sum(mean_n, na.rm = TRUE),
    sum_median_n = sum(median_n, na.rm = TRUE),
    sum_max_n = sum(max_n, na.rm = TRUE),
    sum_median_ncases = sum(median_ncases, na.rm = TRUE),
    sum_median_ncontrol = sum(median_ncontrol, na.rm = TRUE)) %>%
    # round to nearest integers
    mutate_if(is.numeric, round) %>%
    mutate(analysis = paste(study_name, "left out"))
}

###############################################################################
#                            Run MVMR sensitivities                           #
###############################################################################

### parameters
outcomes_list <- model_1_outcomes_list

### maternal GWAS leave-one-study-out & leave out overlapping studies
# read in single study results
stu_res_all <- data.frame()
for(study_name in studies){
    stu_res_tmp_path <- paste0(home_dir, "/results/c_mvmr/",
        study_name, "/raw/mvmr_results_twosample.csv")
    if(file.exists(stu_res_tmp_path)){
        stu_res_tmp <- read.csv(stu_res_tmp_path) %>%
        mutate(analysis = study_name)
        stu_res_all <- rbind(stu_res_all, stu_res_tmp)
    }
}

# calculate sample sizes
loo_sample_sizes <- data.frame()
for(study_name in c(studies, "UKB & ALSPAC")){ 
  loo_tmp <- calculate_loo_sample_size(study_name) %>%
    mutate(exposure = id.exposure)
  loo_sample_sizes <- rbind(loo_sample_sizes, loo_tmp)
}

# analyse each outcome and model exposure in turn
study_loo_res <- data.frame()

# run combined analyses
for(outcome_name in outcomes_list){

    for(exposure_id in c("ieu-b-5136", "ieu-b-5107")){

    # select single study results for exposure-outcome
    stu_res_tmp <- stu_res_all %>%
        filter(id.exposure == exposure_id & outcome == outcome_name)

    ### 1 - leave-one-study out
    study_loo_dat <- metafor::rma(yi = b, sei = se, data = stu_res_tmp,
        slab = analysis, method = "FE")
    # leave one out for all studies
    study_loo_res_tmp <- metafor::leave1out(study_loo_dat) %>%
        as.data.frame %>%
        tibble::rownames_to_column(., var = "study") %>%
        mutate(se = estimate/zval,
            # to remove numbers after study:
            analysis = paste(str_remove(study, "[.][0-9]*"), "left out"),
            exposure = exposure_id,
            outcome = outcome_name)
    study_loo_res <- rbind(study_loo_res, study_loo_res_tmp)

    ### 2 - leave out the overlapping studies specifically
    no_overlap_tmp <- filter(stu_res_tmp, !(analysis %in% overlap_studies))
    # meta-analysis
    no_overlap_res_tmp <- metafor::rma(yi = b, sei = se,
        data = no_overlap_tmp, method = "FE")
    no_overlap_res_tmp_2 <- no_overlap_res_tmp %>%
        tidy() %>% as.data.frame %>%
        # to match first loo analysis columns
        mutate(analysis = "UKB & ALSPAC left out",
            exposure = exposure_id, outcome = outcome_name,
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
    # select same columns as first loo analysis
    no_overlap_res_tmp_2 <- no_overlap_res_tmp_2 %>% 
        select(colnames(study_loo_res_tmp))
    study_loo_res <- rbind(study_loo_res, no_overlap_res_tmp_2)

    }
}

# save results
loo_dir <- paste0(home_dir, "/results/c_mvmr/", analysis, "/loo")
system(paste0("mkdir -p ", loo_dir))
study_loo_res %>%
    left_join(outcome_names_df) %>% 
    left_join(loo_sample_sizes, by = join_by(outcome, analysis, exposure)) %>%
  write.csv(paste0(loo_dir, "/leave_one_study_out_results.csv"), row.names = FALSE)

q('no')