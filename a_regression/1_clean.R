###############################################################################
#                                                                             #
#                       Clean ALSPAC data, generate Table 1                   #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
x <- c("dplyr", "tidyr", "knitr", "data.table", "purrr", "haven",
    "lmtest", "tibble", "skimr", "jtools", "gtsummary")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- Sys.getenv("AGE_AT_MENARCHE_DIR")
setwd(home_dir)
# where dataset is saved
data_dir <- paste0(home_dir, "/working/data/a_regression/")
# to write results to
out_dir <- paste0(home_dir, "/working/results/a_regression/main/")

# Setting digits
options(digits = 10)

# List of pre-selected outcomes
outcome_names_df <- read.csv(file.path("working/data/b_mr/outcome/outcomes.csv"))

###############################################################################
#    Load and filter ALSPAC data on singleton pregnancies to unique mothers   #
###############################################################################

### Load
# De-identified ALSPAC data with all necessary variables
alspac_raw <- read_dta(
    paste0(data_dir, "***.dta"))
# De-identified ALSPAC data with derived MR-PREG outcome variables
df <- read.table(paste0(data_dir, "***.csv"), header = T,
    stringsAsFactors = F, sep = ";")

### Combine these datasets
# Check each pregnancy id and child from raw appears once in each dataset
dim(left_join(df, alspac_raw, by = join_by(cidB4227, qlet)))[1] ==
    dim(inner_join(df, alspac_raw, by = join_by(cidB4227, qlet)))[1]
# Append all ALSPAC variables to cleaned outcome variables
df_extra <- left_join(df, alspac_raw, by = join_by(cidB4227, qlet))

### Filter
dim(df_extra)

# Keep only singleton births
df_extra <- df_extra %>% filter(singleton == 1)
dim(df_extra)

# Keep only unique mothers using in-house variable
df_extra <- df_extra %>% filter(mz005l == 2)
dim(df_extra)

###############################################################################
#                 Construct continuous age at menarche variable               #
###############################################################################

### Clean
# d010
df_extra$d010_clean <- if_else(
    (df_extra$d010 < 0 | df_extra$d010 == 77 | df_extra$d010 == 99), NA,
    df_extra$d010)
table(df_extra$d010, useNA = "always")
table(df_extra$d010_clean, useNA = "always")

# n1120
df_extra$n1120_clean <- if_else(
    (df_extra$n1120 < 0 | df_extra$n1120 == 77 | df_extra$n1120 == 99), NA,
    df_extra$n1120)
table(df_extra$n1120, useNA = "always")
table(df_extra$n1120_clean, useNA = "always")

# r2080
df_extra$r2080_clean <- if_else(
    (df_extra$r2080 < 0 | df_extra$r2080 == 77 | df_extra$r2080 == 99), NA,
    df_extra$r2080)
table(df_extra$r2080, useNA = "always")
table(df_extra$r2080_clean, useNA = "always")

### Combine
# Where participants have responded to this question at multiple questionnaires
# the earliest answer is used
# First, make temporal order explicit
df_extra <- df_extra %>% arrange(d010_clean, n1120_clean, r2080_clean)
# Next, create new variable
df_extra$age_at_menarche <- coalesce(df_extra$d010_clean, df_extra$n1120_clean,
    df_extra$r2080_clean) # (Okay to ignore warnings about value labels here)

### Inspect
table(df_extra$age_at_menarche, useNA = "always")
table(df_extra$age_at_menarche, df_extra$d010_clean)
table(df_extra$age_at_menarche, df_extra$n1120_clean)
table(df_extra$age_at_menarche, df_extra$r2080_clean)

###############################################################################
#                 Clean covariates and baseline characteristics               #
###############################################################################

# Age at delivery, mz028b
df_extra$age_at_del <- if_else((df_extra$mz028b < 0), NA, df_extra$mz028b)

# Offspring sex at birth, kz021
# initially 1 = male, 2 = female; recoded as 0 male & 1 female
df_extra$female_offspring <- if_else((df_extra$kz021 == 2), 1, 0)
# Categoric version as male/female for Table 1
df_extra$offspring_sex <- factor(df_extra$female_offspring,
    labels = c("Male", "Female"))

# Parity, b032
df_extra$parity <- if_else((df_extra$b032 < 0), NA, df_extra$b032)
# All values 4 and over binned; still treated as continuous
df_extra$parity_bins_numeric <- if_else((df_extra$parity >=4), 4, df_extra$parity)
# Categoric version as 0,1,2,3,4+ for Table 1
df_extra$parity_categoric <- factor(df_extra$parity_bins_numeric,
    labels = c("0", "1", "2", "3", "4<="))

# Ever smoked, combining n5000 & r6010
# coded as 1 = ever, 2 = never
df_extra$n5000_clean <- if_else((df_extra$n5000 < 0), NA, df_extra$n5000)
df_extra$r6010_clean <- if_else((df_extra$r6010 < 0), NA, df_extra$r6010)
df_extra <- df_extra %>% arrange(n5000_clean, r6010_clean)
# Categoric version as yes/no/missing for Table 1
df_extra$ever_smoked <-
    coalesce(df_extra$n5000_clean, df_extra$r6010_clean)
df_extra$ever_smoked[is.na(df_extra$ever_smoked)] <- 100
df_extra$ever_smoked <- factor(df_extra$ever_smoked,
    labels = c("Yes", "No", "Missing"))

# Alcohol consumption before this pregnancy, b720
df_extra <- df_extra %>%
    mutate(alcohol = case_when(b720 < 0 ~ NA, b720 == 9 ~ NA, b720 == 6 ~ 5,
        .default = b720))
df_extra$alcohol[is.na(df_extra$alcohol)] <- 100
# Categoric version incluing missing for Table 1
df_extra$alcohol <- factor(df_extra$alcohol,
    labels = c("Never", "<1 glass per week", "1+ glass per week",
        "1-2 glasses per day", "3+ glasses per day", "Missing"))

# Mother's highest educational qualification (categoric, 5), c645a
attr(df_extra$c645a, "labels")
df_extra$highest_ed <- if_else((df_extra$c645a < 0), NA, df_extra$c645a)
df_extra$highest_ed <- factor(df_extra$highest_ed,
  labels = c("CSE", "Vocational", "O-level", "A-level", "Degree"))

# Ethnicity, code as 1 = white & 0 = non-white, c800
df_extra$ethnicity <- if_else(
  (df_extra$c800 < 0 | df_extra$c800 == 99), NA, df_extra$c800)
df_extra$ethnicity <- if_else((df_extra$ethnicity == 1), 1, 0)
df_extra$ethnicity <- factor(df_extra$ethnicity,
 labels = c("non-White", "White"))

# BMI at 12 weeks' gestation, dw042
df_extra$bmi <- if_else((df_extra$dw042 < 0), NA, df_extra$dw042)

###############################################################################
#                           Filter to complete cases                          #
###############################################################################

dim(df_extra)

### Filter out participants missing age at menarche
df_exposure_complete <- df_extra %>% drop_na(age_at_menarche)
dim(df_exposure_complete)

### Filter out participants missing all outcomes
# (i.e. keep all with data on at least one)
df_exp_out_complete <- df_exposure_complete %>%
    filter_at(vars(any_of(all_outcomes)), any_vars(!(is.na(.))))
dim(df_exp_out_complete)

### Filter out participants missing any model covariates
# First, specify model covariates
covariates_vector <- c("age_at_del", "parity_bins_numeric", "female_offspring",
    "highest_ed", "ethnicity")
covariates_bmi_vector <- c(covariates_vector, "bmi")
# Second, filter
df_complete <- df_exp_out_complete %>% drop_na(all_of(covariates_bmi_vector))
dim(df_complete)

###############################################################################
#                  Construct age at menarche variables                        #
###############################################################################

### Construct early, normative and late categories based on the sample mean and
# SD
mean_aam <- mean(df_complete$age_at_menarche, na.rm = TRUE)
print(mean_aam)
sd_aam <- sd(df_complete$age_at_menarche, na.rm = TRUE)
print(sd_aam)

# Construct
df_complete <- df_complete %>%
  mutate(age_at_menarche_cat =
  case_when((age_at_menarche > (mean_aam - sd_aam) &
                age_at_menarche < (mean_aam + sd_aam)) ~ "normative", # 12 to 14
            age_at_menarche < (mean_aam - sd_aam) ~ "early", # 11 and under
            age_at_menarche > (mean_aam + sd_aam) ~ "late", # 15 and over
            .default = NA))
# Order as factor
df_complete$age_at_menarche_cat <- factor(df_complete$age_at_menarche_cat,
    ordered = TRUE, levels = c("normative", "early", "late"))
# Variation of the variable for modelling as numeric 1/2/3 for LRTest
df_complete <- df_complete %>%
  mutate(age_at_menarche_cat_numeric =
  case_when(age_at_menarche_cat == "early" ~ 1, # 11 and under
            age_at_menarche_cat == "normative" ~ 2,
            age_at_menarche_cat == "late" ~ 3, # 15 and over
            .default = NA))

### Bin extreme values in continuous version
# Group 9 and under (8,9) and 17 and over (17,18,...)
df_complete <- df_complete %>%
    mutate(age_at_menarche_bins = cut(age_at_menarche,
    breaks = c(-Inf,9:16,Inf)))
# Inspect
table(df_complete$age_at_menarche_bins, df_complete$age_at_menarche)
# Continuous version of variable
df_complete <- df_complete %>%
  mutate(age_at_menarche_bins_numeric = as.numeric(age_at_menarche_bins))

###############################################################################
#                                      Table 1                                #
###############################################################################

# Describe whole sample:
# age at delivery
mean(df_complete$age_at_del); sd(df_complete$age_at_del)
# age at menarche
mean(df_complete$age_at_menarche); sd(df_complete$age_at_menarche)

# Table 1 using gtsummary package, by  1 = early, 2 = normative, 3 = late
table1 <- df_complete %>% 
  select(age_at_menarche_cat_numeric, age_at_del, bmi, offspring_sex,
    parity_categoric, ever_smoked, alcohol, highest_ed, ethnicity) %>%
  tbl_summary(
    by = age_at_menarche_cat_numeric,
    statistic = list(all_categorical() ~ "{n} ({p}%)",
                     age_at_del ~ "{mean} ({sd})",
                     bmi ~ "{mean} ({sd})"),
    digits = list(all_continuous() ~ c(2, 2),
                  all_categorical() ~ c(0, 1)),
    type = list(age_at_menarche_cat_numeric ~ "categorical",
                age_at_del ~ "continuous",
                bmi ~ "continuous",
                offspring_sex ~ "categorical",
                parity_categoric ~ "categorical",
                ever_smoked ~ "categorical",
                alcohol ~ "categorical",
                highest_ed ~ "categorical",
                ethnicity ~ "categorical"),
    label = list(age_at_menarche_cat_numeric ~ "Menarche onset",
                age_at_del ~ "Maternal age at delivery",
                bmi ~ "BMI at 12 weeks' gestation",
                offspring_sex ~ "Offspring sex",
                parity_categoric ~ "Parity",
                ever_smoked ~ "Ever smoked",
                alcohol ~ "Alcohol before pregnancy",
                highest_ed ~ "Education",
                ethnicity ~ "Ethnicity")
            ) %>%
    as.data.frame()

print(table1)
write.csv(table1,
    paste0(out_dir, "table_1.csv"),
    row.names = FALSE)

###############################################################################
#                                Save cleaned data                            #
###############################################################################

# variables needed for modelling
final_variables <- c(all_outcomes, "age_at_menarche", "age_at_del",
    "offspring_sex", "parity_bins_numeric", "highest_ed", "ethnicity",
    "bmi", "age_at_menarche_cat", "age_at_menarche_cat_numeric",
    "age_at_menarche_bins", "age_at_menarche_bins_numeric")

df_final <- df_complete %>%
    select(all_of(final_variables))

write.csv(df_final, paste0(data_dir, "/cleaned/df_final.csv"),
    row.names = FALSE)

q('no')