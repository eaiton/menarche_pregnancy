###############################################################################
#                                                                             #
#                     Calculate sample overlap in GWAS                        #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear environment
rm(list = ls())

# Required libraries
x <- c("dplyr", "tidyr", "data.table")
lapply(x, require, character.only = TRUE)

library(TwoSampleMR)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "working/")
out_dir <- file.path(home_dir, "results/d_combined/sample_overlap/")
data_dir <- paste0(home_dir, "data/b_mr/")
outcome_data_dir <- file.path(data_dir, "outcome/")
setwd(home_dir)

univariate_mr_binary_path <- paste0(home_dir, "results/b_mr/mr_results_binary.csv")
univariate_mr_bw_path <- paste0(home_dir, "results/b_mr/mr_results_bw.csv")
mvmr_binary_path <- paste0(home_dir, "results/c_mvmr/mvmr_results_binary.csv")
mvmr_bw_path <- paste0(home_dir, "results/c_mvmr/mvmr_results_bw.csv")

# Outcome names
outcome_names_df <- read.csv(file.path(outcome_data_dir, "/outcomes.csv"))
outcomes <- c(outcome_names_df$outcome)

###############################################################################
#                        Exposures GWAS sample sizes                          #
###############################################################################

# Age at menarche - Kentistou et al 2017 - total size 632,955
kentistou_aam_total <- 632955
# Study-level sample sizes from supplementary table S1
kentistou_aam_alspac <- 9283
# ALSPAC both mothers and children were used
# but it is not clear the exact numbers from each generation
# But only the mothers generation (G0) would overlap with outcome GWAS
# Conservatively assuming 100% are mothers for ease of calculation
kentistou_aam_ukbb <- 238040
kentistou_aam_overlap <- kentistou_aam_alspac + kentistou_aam_ukbb

kentistou_aam_overlap/kentistou_aam_total*100
# 39% of exposure GWAS is potentially in outcome GWAS

# Childhood body size - Richardson et al 2020 - total size 453,169
richardson_adiposity_total <- 453169
# All UKB, but only females in sample will overlap pregnancy outcome data
richardson_adiposity_overlap <- 246511

###############################################################################
#                   Outcome GWAS sample sizes - by outcome                    #
###############################################################################

### Univariate MR analyses total meta-analysis sample sizes
### using medians since totals vary by SNP
mr_binary <- read.csv(univariate_mr_binary_path) %>%
  filter(method == "Inverse variance weighted", analysis == "Main") %>%
  mutate(outcome_gwas_total = median_n) %>% 
  select(outcome_clean, outcome_gwas_total)
mr_bw <- read.csv(univariate_mr_bw_path) %>% 
  filter(method == "Inverse variance weighted", analysis == "Main") %>%
  mutate(outcome_gwas_total = median_n) %>% 
  select(outcome_clean, outcome_gwas_total)
mr <- rbind(mr_binary, mr_bw) %>% 
  arrange(desc(outcome_gwas_total))

# MVMR analyses total meta-analysis sample sizes
### using medians since totals vary by SNP
mvmr_binary <- read.csv(mvmr_binary_path) %>%
  filter(analysis == "Main", id.exposure == "ieu-b-5136") %>% 
  mutate(outcome_gwas_total = median_n) %>% 
  select(outcome_clean, outcome_gwas_total)
mvmr_bw <- read.csv(mvmr_bw_path) %>% 
  filter(analysis == "Main", id.exposure == "ieu-b-5136") %>% 
  mutate(outcome_gwas_total = median_n) %>% 
  select(outcome_clean, outcome_gwas_total)
mvmr <- rbind(mvmr_binary, mvmr_bw) %>% 
  arrange(desc(outcome_gwas_total))

# mr and mvmr are almost identical numbers
# differences due to different SNPs used

###############################################################################
#             Exposure GWAS potential overlap - by outcome                    #
###############################################################################

### AAM GWAS with outcome GWAS
### Applies to univariate and MVMR
### Median size of UKB & ALSPAC GWAS for each outcome
mr_binary_overlap <- read.csv(univariate_mr_binary_path) %>%
  filter(analysis == "ALSPAC only" | analysis == "UKB only") %>%
  filter(method == "Inverse variance weighted") %>% 
  group_by(outcome_clean) %>%
  summarise(overlap_gwas_total = sum(median_n)) %>%
  select(outcome_clean, overlap_gwas_total)
mr_bw_overlap <- read.csv(univariate_mr_bw_path) %>% 
  filter(analysis == "ALSPAC only" | analysis == "UKB only") %>%
  filter(method == "Inverse variance weighted") %>% 
  group_by(outcome_clean) %>%
  summarise(overlap_gwas_total = sum(median_n)) %>%
  select(outcome_clean, overlap_gwas_total)
mr_overlap <- rbind(mr_binary_overlap, mr_bw_overlap) %>% 
 arrange(desc(overlap_gwas_total)) 

mr_overlap$overlap_gwas_total < kentistou_aam_total
# since this is true for all, can use this as maximal possible overlap
# if false, would use kentistou_aam_total instead

### Adiposity GWAS with outcome GWAS
### MVMR only
### Median size of UKB females for each outcome
### combined male & female participants used in exposure GWAS,
### but outcome GWAS restricted to female participants
female_prop <- richardson_adiposity_overlap/richardson_adiposity_total
# ~54%
mvmr_binary_overlap_adiposity <- read.csv(mvmr_binary_path) %>%
  filter(analysis == "UKB only") %>%
  filter(id.exposure == "ieu-b-5107") %>% # doesn't matter
  mutate(median_n_female = median_n * female_prop) %>%
  group_by(outcome_clean) %>%
  summarise(overlap_gwas_adiposity = sum(median_n_female)) %>%
  select(outcome_clean, overlap_gwas_adiposity)
mvmr_bw_overlap_adiposity <- read.csv(mvmr_bw_path) %>%
  filter(analysis == "UKB only") %>%
  filter(id.exposure == "ieu-b-5107") %>% # doesn't matter
  mutate(median_n_female = median_n * female_prop) %>%
  group_by(outcome_clean) %>%
  summarise(overlap_gwas_adiposity = sum(median_n_female)) %>%
  select(outcome_clean, overlap_gwas_adiposity)
mvmr_overlap_adiposity <- rbind(mvmr_binary_overlap_adiposity, mvmr_bw_overlap_adiposity) %>% 
 arrange(desc(overlap_gwas_adiposity))
# this is maximal possible overlap

###############################################################################
#                      Calculate overlap - by outcome                         #
###############################################################################

### age at menarche maximum overlap per outcome - univariate MR
aam_mr_overlap <- mr %>% 
    left_join(mr_overlap) %>%
    mutate(exposure_gwas_total = kentistou_aam_total,
          # maximum possible overlap -
           exposure_outcome_overlap = overlap_gwas_total,
           percent_exposure_in_outcome_gwas = exposure_outcome_overlap/outcome_gwas_total*100,
           percent_outcome_in_exposure_gwas = exposure_outcome_overlap/exposure_gwas_total*100)
aam_mr_overlap
write.csv(aam_mr_overlap,
  paste0(out_dir, "aam_mr_overlap.csv"), quote = T, row.names = F)

### age at menarche maximum overlap per outcome - MVMR
aam_mvmr_overlap <- mvmr %>% 
    left_join(mr_overlap) %>%
    mutate(exposure_gwas_total = kentistou_aam_total,
          # maximum possible overlap -
           exposure_outcome_overlap = overlap_gwas_total,
           percent_exposure_in_outcome_gwas = exposure_outcome_overlap/outcome_gwas_total*100,
           percent_outcome_in_exposure_gwas = exposure_outcome_overlap/exposure_gwas_total*100)
aam_mvmr_overlap
write.csv(aam_mvmr_overlap,
  paste0(out_dir, "aam_mvmr_overlap.csv"), quote = T, row.names = F)

### adiposity aged 10 maximum overlap per outcome - MVMR
adiposity_mvmr_overlap <- mvmr %>% 
    left_join(mvmr_overlap_adiposity) %>%
    mutate(exposure_gwas_total = richardson_adiposity_total,
          # maximum possible overlap -
           exposure_outcome_overlap = overlap_gwas_adiposity,
           percent_exposure_in_outcome_gwas = exposure_outcome_overlap/outcome_gwas_total*100,
           percent_outcome_in_exposure_gwas = exposure_outcome_overlap/exposure_gwas_total*100)
adiposity_mvmr_overlap
write.csv(adiposity_mvmr_overlap,
  paste0(out_dir, "adiposity_mvmr_overlap.csv"), quote = T, row.names = F)

q('no')