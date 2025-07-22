###############################################################################
#                                                                             #
#                Download GWAS summary statistics locally                     #
#                      for all MR & MVMR exposures                            #
#                                                                             #
###############################################################################

###############################################################################
#                                  Set up                                     #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
# remotes::install_github('MRCIEU/TwoSampleMR')
# remotes::install_github("MRCIEU/MRInstruments")
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "working/")
rdsf_ieu2 <- Sys.getenv("RDSF_IEU2")
# Where to save GWAS VCF files
out_dir <- paste0(home_dir, "data/b_mr/gwas/")
setwd(home_dir)

###############################################################################
#                              List exposures                                 #
###############################################################################

### all univariate & multivariate MR exposures

## make csv
phenotype <- c("Age at menarche (years)",
    "Comparative body size aged 10 (SD)")
phenotype_id <- c("aam", "body-size-10")
opengwas_id <- c("ieu-b-5136", "ieu-b-5107")
sample_size <- c(632955, 453169)
opengwas_list <- data.frame(cbind(phenotype, phenotype_id, opengwas_id, sample_size))
write.csv(opengwas_list, paste0(out_dir, "/opengwas_list.csv"),
    row.names = FALSE)

###############################################################################
#              Check if GWAS are on IEU Open GWAS using API                   #
###############################################################################

instr <- extract_instruments(outcomes = opengwas_id)
table(instr$exposure)

###############################################################################
#                   Download VCFs for all traits locally                      #
###############################################################################

# so that GWAS summary statistics can then be accessed without relying on
# IEU OpenGWAS API servers, see: https://gwas.mrcieu.ac.uk/about/

# use csv of traits & codes 
for(id in opengwas_id){
    file_url = paste0("https://gwas.mrcieu.ac.uk/files/", id, "/", id, ".vcf.gz")
    system(paste0("wget -nc ", file_url, " -O ", out_dir, id, ".vcf.gz"))
}

q('no')