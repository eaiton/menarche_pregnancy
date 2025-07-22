###############################################################################
#                        Select exposure SNPs for MVMR                        #
#                           using TwoSampleMR API                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list=ls())

# Required libraries
# remotes::install_github('MRCIEU/TwoSampleMR')
library(TwoSampleMR)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "/working/")
gwas_dir <- file.path(home_dir, "/data/b_mr/gwas/")
out_dir <- file.path(home_dir, "/data/c_mvmr/")

# List of exposures
opengwas_list <- read.csv(paste0(gwas_dir, "/opengwas_list.csv"))
all_opengwasid <- c(opengwas_list$opengwas_id)

# Check token
ieugwasr::get_opengwas_jwt()
Sys.setenv(OPENGWAS_JWT="***")

###############################################################################
#                           Specify model exposures                           #
###############################################################################

### Model 1 - age at menarche & body size aged 10
model_1_ids <- all_opengwasid[c(1,3)]

###############################################################################
#                         Retrieve clumped instruments                        #
###############################################################################

# Using TwoSampleMR function to query API
# Default function settings match local clumping used in univariate MR:
#  clump_r2 = 0.001, clump_kb = 10000, find_proxies = TRUE

exposure_dat_model_1 <- mv_extract_exposures(model_1_ids)
dim(exposure_dat_model_1)
table(exposure_dat_model_1$exposure)

exposure_dat_model_1 <- mv_extract_exposures(model_1_ids[2], model_1_ids[1])
dim(exposure_dat_model_1)
table(exposure_dat_model_1$exposure)

###############################################################################
#                             Instrument details                              #
###############################################################################

### Instruments for each exposure
aam <- extract_instruments("ieu-b-5136")
dim(aam) # 467
body_size <- extract_instruments("ieu-b-5107")
dim(body_size) # 231

### Unique instruments
snps <- c(aam$SNP, body_size$SNP)
length(snps) # 698
length(unique(snps)) # 696

###############################################################################
#                  Write out clumped exposure instruments                     #
###############################################################################

write.csv(exposure_dat_model_1,
  paste0(out_dir, "instruments/mvmr_snps.csv"),
  quote = T, row.names = F)

q('no')