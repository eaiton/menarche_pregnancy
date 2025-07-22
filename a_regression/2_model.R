###############################################################################
#                                                                             #
#                     Multivariable regression in ALSPAC                      #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
x <- c("dplyr", "tidyr", "knitr", "data.table", "purrr", "lmtest", "tibble",
  "jtools", "gtsummary", "broom")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- Sys.getenv("AGE_AT_MENARCHE_DIR")
setwd(home_dir)
# where dataset is saved
data_dir <- paste0(home_dir, "/working/data/a_regression/")
# to write results to
out_dir <- paste0(home_dir, "/working/results/a_regression/main/")

# Ensure models are run using complete cases
options(na.action = "na.omit")

# Setting digits
options(digits = 10)

# Outcomes
outcome_names_df <- read.csv(file.path("working/data/b_mr/outcome/outcomes.csv"))
all_outcomes <- outcome_names_df$outcome
binary_outcomes <- all_outcomes[all_outcomes != "zbw_all"]
continuous_outcomes <- c("zbw_all")

# Set model covariates
# all covariates
covariates_vector <- c("age_at_del", "parity_bins_numeric", "offspring_sex",
  "highest_ed", "ethnicity")
covariates <-
  " + age_at_del + parity_bins_numeric + offspring_sex + highest_ed + ethnicity"
# all covariates and bmi
covariates_bmi_vector <- c(covariates_vector, "bmi")
covariates_bmi <- paste0(covariates, " + bmi")

###############################################################################
#                           Load cleaned ALSPAC data                          #
###############################################################################

### Load
df <- read.csv(paste0(data_dir, "/cleaned/df_final.csv"))

### Ensure variable types as expected for modelling
# Recode all binary categoric variables as factors
factors <- c("offspring_sex", "ethnicity")
df[, factors] <- lapply(df[, factors] , factor)

# Recode categoric variables with several levels in correct order
# Highest education
df$highest_ed_factor <- factor(df$highest_ed,
  levels = c("CSE", "Vocational", "O-level", "A-level", "Degree"))
# Age at menarche categories
df$age_at_menarche_cat <- factor(df$age_at_menarche_cat,
  levels = c("early", "normative", "late"))
df$age_at_menarche_bins <- factor(df$age_at_menarche_bins,
    levels = c("(-Inf,9]", "(10,11]", "(11,12]", "(12,13]", "(13,14]",
    "(14,15]", "(15,16]", "(16, Inf]"))
# Generate two factor variable for age at menarche categories:
# early (1) as baseline for LRT tests
df$age_at_menarche_cat_numeric_factor <- factor(df$age_at_menarche_cat_numeric,
    levels = c(1:3))
# Re-levelled with normative (2) as baseline for categoric models
df$age_at_menarche_cat_relevel <-
  relevel(as.factor(df$age_at_menarche_cat_numeric), ref = 2)

### Inspect
str(df)

###############################################################################
#                      LRT of linearity for each outcome                      #
###############################################################################

# Comparing nested models with age at menarche exposure as early-normative-late
# encoded as 1/2/3 categories compared to 1->2->3 continuous

lrt_results <- data.frame()
lrt_results_raw <- data.frame()

for(i in seq_along(all_outcomes)) {

    outcome <- all_outcomes[i]
  
    ### Fit nested models with age at menarche as numeric and as categoric
    if(outcome %in% binary_outcomes){
      tmp_linear_model <-
          glm(paste0(outcome, " ~ age_at_menarche_cat_numeric", covariates_bmi),
              df, family = "binomial")
      tmp_categoric_model <-
          glm(paste0(outcome, " ~ age_at_menarche_cat_numeric_factor",
              covariates_bmi), df, family = "binomial")
      } else {
      tmp_linear_model <-
          glm(paste0(outcome, " ~ age_at_menarche_cat_numeric", covariates_bmi),
              df, family = "gaussian")
      tmp_categoric_model <-
          glm(paste0(outcome, " ~ age_at_menarche_cat_numeric_factor",
              covariates_bmi), df, family = "gaussian")
    }

    ### Likelihood ratio test to assess goodness-of-fit
    tmp_lrt <- lrtest(tmp_linear_model, tmp_categoric_model)
    tmp_lrt_tidy <- tidy(tmp_lrt) %>% 
      mutate(LogLik_rounded = round(LogLik, 2),
             Chi_rounded = round(statistic, 2),
             model = ifelse(grepl("age_at_menarche_cat_numeric_factor",
                                        term, ignore.case = FALSE), "Categoric", "Linear"),
             outcome = outcome)
    lrt_results_raw <- rbind(lrt_results_raw, tidy(tmp_lrt))
    
    lrt_results <- rbind(lrt_results, tmp_lrt_tidy)

}

print(lrt_results)
print(lrt_results_raw)

lrt_results %>% 
  left_join(outcome_names_df) %>% 
  select(outcome_clean, model, LogLik_rounded, X.Df, Chi_rounded, df, p.value) %>% 
  pivot_wider(names_from = model, values_from = c(LogLik_rounded, X.Df)) %>% 
  write.csv(paste0(out_dir, "/lrt_results.csv"),
    row.names = FALSE)

write.csv(lrt_results_raw, paste0(out_dir, "/lrt_results_raw.csv"),
          row.names = FALSE)

###############################################################################
#                  Functions to fit and summarise models                      #
###############################################################################

### Define functions to fit all models and summarise, allowing for both binary
### & continuous exposures

### Binary outcomes function
fit_and_summarise_binomial_model <-
  function(model_specified, df, outcome, exposure_variable, model_name) {

  ### Filter to complete cases for this outcome
  model_cols <- c(exposure_variable, covariates_bmi_vector, outcome)
  model_complete_cases <- drop_na(df, all_of(model_cols))
  
  ### Fit appropriate regression model for outcome
  ### Check using a binary outcome, OR
  if(outcome %in% binary_outcomes){
  model <- glm(model_specified, model_complete_cases, family = "binomial")
  ### b. Continuous outcome, beta
  } else {
    warning("Outcome not in binary outcome list")
  }
  
  ### Summarise model
  ### a. Age at menarche as a linear exposure
  if(exposure_variable == "age_at_menarche_bins_numeric"){
  model_out <- model %>%
    summ(confint = TRUE, digits = 10, exp = TRUE) %>%
    tidy() %>%
    mutate(outcome = outcome,
           model_specified = model_specified,
           model_name = model_name,
           # OR scale
           beta_OR = estimate,
           std.error_OR = std.error,
           # log OR scale
           beta_logOR = log(beta_OR),
           std.error_logOR = log(std.error_OR))
  # Sample size
  model_out$model_nobs <- nobs(model)
  model_out$model_ncases <- sum(model_complete_cases[, outcome])
  model_out$model_ncontrols <- model_out$model_nobs - model_out$model_ncases
  # Select only exposure coefficient
  model_out <- filter(model_out, term == exposure_variable)

  ### b. Age at menarche as a categoric exposure
  } else if (exposure_variable == "age_at_menarche_cat_relevel") {
  model_out <- model %>%
    summ(confint = TRUE, digits = 10, exp = TRUE) %>%
    tidy() %>%
    mutate(outcome = outcome,
           model_specified = model_specified,
           model_name = model_name)
  # Estimates saved as new columns, wide format
  # logOR scale
  model_out$beta_logOR_norm_vs_early <-
    summary(model)$coefficients["age_at_menarche_cat_relevel1", "Estimate"]
  model_out$beta_logOR_norm_vs_late <-
    summary(model)$coefficients["age_at_menarche_cat_relevel3", "Estimate"]
  model_out$std.error_logOR_norm_vs_early <-
    summary(model)$coefficients["age_at_menarche_cat_relevel1", "Std. Error"]
  model_out$std.error_logOR_norm_vs_late <-
    summary(model)$coefficients["age_at_menarche_cat_relevel3", "Std. Error"]
  # OR scale, just exp() above ones
  model_out$beta_OR_norm_vs_early <-
    exp(model_out$beta_logOR_norm_vs_early)
  model_out$beta_OR_norm_vs_late <-
    exp(model_out$beta_logOR_norm_vs_late)
  model_out$std.error_OR_norm_vs_early <-
    exp(model_out$std.error_logOR_norm_vs_early)
  model_out$std.error_OR_norm_vs_late <-
    exp(model_out$std.error_logOR_norm_vs_late)
  # Sample size
  model_out$model_nobs <- nobs(model)
  model_out$model_ncases <- sum(model_complete_cases[, outcome])
  model_out$model_ncontrols <- model_out$model_nobs - model_out$model_ncases
  # Select only exposure coefficients
  model_out <- filter(model_out, (term == "age_at_menarche_cat_relevel1"
    | term == "age_at_menarche_cat_relevel3"))
  }
  return(model_out)
}

### Continuous outcomes function
fit_and_summarise_linear_model <-
  function(model_specified, df, outcome, exposure_variable, model_name) {

  ### Filter to complete cases for this outcome
  model_cols <- c(exposure_variable, covariates_bmi_vector, outcome)
  model_complete_cases <- drop_na(df, all_of(model_cols))
  
  ### Fit appropriate regression model for outcome
  ### Check using a binary outcome, OR
  if(outcome %in% continuous_outcomes){
  model <- glm(model_specified, model_complete_cases, family = "gaussian")
  ### b. Continuous outcome, beta
  } else {
    warning("Outcome not in continuous outcome list")
  }
  
  ### Summarise model
  ### a. Age at menarche as a linear exposure
  if(exposure_variable == "age_at_menarche_bins_numeric"){
  model_out <- model %>%
    summ(confint = TRUE, digits = 10, exp = FALSE) %>%
    tidy() %>%
    mutate(outcome = outcome,
           model_specified = model_specified,
           model_name = model_name,
           # betas
           beta = estimate,
           std.error = std.error)
  # Sample size
  model_out$model_nobs <- nobs(model)
  # Select only exposure coefficient
  model_out <- filter(model_out, term == exposure_variable)

  ### b. Age at menarche as a categoric exposure
  } else if (exposure_variable == "age_at_menarche_cat_relevel") {
  model_out <- model %>%
    summ(confint = TRUE, digits = 10, exp = FALSE) %>%
    tidy() %>%
    mutate(outcome = outcome,
           model_specified = model_specified,
           model_name = model_name)
  # Estimates saved as new columns, wide format
  # beta scale
  model_out$beta_norm_vs_early <-
    summary(model)$coefficients["age_at_menarche_cat_relevel1", "Estimate"]
  model_out$beta_norm_vs_late <-
    summary(model)$coefficients["age_at_menarche_cat_relevel3", "Estimate"]
  model_out$std.error_norm_vs_early <-
    summary(model)$coefficients["age_at_menarche_cat_relevel1", "Std. Error"]
  model_out$std.error_norm_vs_late <-
    summary(model)$coefficients["age_at_menarche_cat_relevel3", "Std. Error"]
  # Sample size
  model_out$model_nobs <- nobs(model)
  # Select only exposure coefficients
  model_out <- filter(model_out, (term == "age_at_menarche_cat_relevel1" |
    term == "age_at_menarche_cat_relevel3"))
  }
  return(model_out)
}

###############################################################################
#                  Primary analysis - linear exposure models                  #
###############################################################################

# Run unadjusted, adjusted, adjusted + adiposity models with age at menarche
# as a continuous exposure (1 unit = 1 year older), as primary analysis

### Binary outcomes

# Create an empty list to store model summaries
linear_exp_binary_out <- list()

# Iterate over binary outcomes
for (outcome in binary_outcomes) {

  # Specify models
  model_raw <- paste0(outcome, " ~ age_at_menarche_bins_numeric")
  model_adjusted <- paste0(outcome, " ~ age_at_menarche_bins_numeric",
    covariates)
  model_adjusted_bmi <- paste0(outcome, " ~ age_at_menarche_bins_numeric",
    covariates_bmi)

  # Fit and summarise binomial models for each
  model_raw_summary <- fit_and_summarise_binomial_model(
    model_raw, df, outcome, "age_at_menarche_bins_numeric", "raw")
  model_adjusted_summary <- fit_and_summarise_binomial_model(
    model_adjusted, df, outcome, "age_at_menarche_bins_numeric", "adjusted")
  model_adjusted_bmi_summary <- fit_and_summarise_binomial_model(
    model_adjusted_bmi, df, outcome, "age_at_menarche_bins_numeric",
    "adjusted_bmi")
  
  # Append model summaries to the list
  linear_exp_binary_out <- c(linear_exp_binary_out, list(model_raw_summary,
    model_adjusted_summary, model_adjusted_bmi_summary))
}

# Combine all model summaries into a single data frame
linear_exp_binary_out <- do.call(rbind, linear_exp_binary_out) %>%
    left_join(outcome_names_df)
write.csv(linear_exp_binary_out, paste0(out_dir, "/linear_exp_binary_out.csv"),
    row.names = FALSE)

### Continuous outcomes

# Create an empty list to store model summaries
linear_exp_cont_out <- list()

# Iterate over continuous outcomes
for (outcome in continuous_outcomes) {

  # Specify models
  model_raw <- paste0(outcome, " ~ age_at_menarche_bins_numeric")
  model_adjusted <- paste0(outcome, " ~ age_at_menarche_bins_numeric",
    covariates)
  model_adjusted_bmi <- paste0(outcome, " ~ age_at_menarche_bins_numeric",
    covariates_bmi)

  # Fit and summarise linear models for each
  model_raw_summary <- fit_and_summarise_linear_model(
    model_raw, df, outcome, "age_at_menarche_bins_numeric", "raw")
  model_adjusted_summary <- fit_and_summarise_linear_model(
    model_adjusted, df, outcome, "age_at_menarche_bins_numeric", "adjusted")
  model_adjusted_bmi_summary <- fit_and_summarise_linear_model(
    model_adjusted_bmi, df, outcome, "age_at_menarche_bins_numeric",
    "adjusted_bmi")
  
  # Append model summaries to the list
  linear_exp_cont_out <- c(linear_exp_cont_out,
    list(model_raw_summary, model_adjusted_summary, model_adjusted_bmi_summary))
}

# Combine all model summaries into a single data frame
linear_exp_cont_out <- do.call(rbind, linear_exp_cont_out) %>%
    left_join(outcome_names_df)
write.csv(linear_exp_cont_out, paste0(out_dir, "/linear_exp_cont_out.csv"),
    row.names = FALSE)

###############################################################################
#                    Sensitivity - categoric exposure models                  #
###############################################################################

# Run unadjusted, adjusted, adjusted + adiposity models with age at menarche
# modelled as normative vs early, normative vs late categories as sensitivity

### Binary outcomes
# Create an empty list to store model summaries
cat_exp_binary_out <- list()

# Iterate over binary outcomes
for (outcome in binary_outcomes) {

  # Specify models
  model_raw <- paste0(outcome, " ~ age_at_menarche_cat_relevel")
  model_adjusted <- paste0(outcome, " ~ age_at_menarche_cat_relevel",
    covariates)
  model_adjusted_bmi <- paste0(outcome, " ~ age_at_menarche_cat_relevel",
    covariates_bmi)

  # Fit and summarise binomial models for each
  model_raw_summary <- fit_and_summarise_binomial_model(
    model_raw, df, outcome, "age_at_menarche_cat_relevel", "raw")
  model_adjusted_summary <- fit_and_summarise_binomial_model(
    model_adjusted, df, outcome, "age_at_menarche_cat_relevel", "adjusted")
  model_adjusted_bmi_summary <- fit_and_summarise_binomial_model(
    model_adjusted_bmi, df, outcome, "age_at_menarche_cat_relevel",
    "adjusted_bmi")
  
  # Append model summaries to the list
  cat_exp_binary_out <- c(cat_exp_binary_out,
    list(model_raw_summary, model_adjusted_summary, model_adjusted_bmi_summary))
}

# Combine all model summaries into a single data frame
cat_exp_binary_out <- do.call(rbind, cat_exp_binary_out) %>%
    left_join(outcome_names_df)
write.csv(cat_exp_binary_out, paste0(out_dir, "/cat_exp_binary_out.csv"),
    row.names = FALSE)

### Continous outcomes
# Create an empty list to store model summaries
cat_exp_cont_out <- list()

# Iterate over continuous outcomes
for (outcome in continuous_outcomes) {

  # Specify models
  model_raw <- paste0(outcome, " ~ age_at_menarche_cat_relevel")
  model_adjusted <- paste0(outcome, " ~ age_at_menarche_cat_relevel",
    covariates)
  model_adjusted_bmi <- paste0(outcome, " ~ age_at_menarche_cat_relevel",
    covariates_bmi)

  # Fit and summarise linear models for each
  model_raw_summary <- fit_and_summarise_linear_model(
    model_raw, df, outcome, "age_at_menarche_cat_relevel", "raw")
  model_adjusted_summary <- fit_and_summarise_linear_model(
    model_adjusted, df, outcome, "age_at_menarche_cat_relevel", "adjusted")
  model_adjusted_bmi_summary <- fit_and_summarise_linear_model(
    model_adjusted_bmi, df, outcome, "age_at_menarche_cat_relevel",
    "adjusted_bmi")
  
  # Append model summaries to the list
  cat_exp_cont_out <- c(cat_exp_cont_out,
    list(model_raw_summary, model_adjusted_summary, model_adjusted_bmi_summary))
}

# Combine all model summaries into a single data frame
cat_exp_cont_out <- do.call(rbind, cat_exp_cont_out) %>%
    left_join(outcome_names_df)
write.csv(cat_exp_cont_out, paste0(out_dir, "/cat_exp_cont_out.csv"),
    row.names = FALSE)

q('no')