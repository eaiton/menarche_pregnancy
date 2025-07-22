###############################################################################
#                                                                             #
#             Save list of all SNPs to extract SNP-outcome data for           #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list = ls())

# Required libraries
x <- c("dplyr", "tibble")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "working/")
gwas_dir <- file.path(home_dir, "data/b_mr/gwas")
out_dir <- file.path(home_dir, "data/b_mr/exposure_dat")
proxies_dir <- file.path(home_dir, "data/_mr/proxies")
setwd(home_dir)

###############################################################################
#                     Read in all instruments and proxies                     #
###############################################################################

all_snps <- c()
all_proxies <- c()

# Read in exposure GWAS information
opengwas_list <- read.csv(file.path(gwas_dir, "opengwas_list.csv"))

for (i in seq_along(opengwas_list$phenotype_id)) {
    phenotype_id <- opengwas_list$phenotype_id[i]
    opengwas_id <- opengwas_list$opengwas_id[i]

    print(paste0("Phenotype: ", phenotype_id, ", ID: ", opengwas_id))

    # Read in clumped genome-wide significant SNPs
    clumped_path <- file.path(out_dir, paste0(phenotype_id, "_pval_5e_08_clumped.csv"))
    tmp_clumped <- read.csv(clumped_path)
    head(tmp_clumped)

    # Read in proxies, r2 > 0.8
    proxies_path <- file.path(out_dir, paste0(phenotype_id, "_proxies.csv"))
    proxies <- read.csv(proxies_path)
    head(proxies)

    all_proxies <- rbind(all_proxies, proxies)

    all_snps <- c(all_snps, tmp_clumped$SNP, proxies$rsid)
}

###############################################################################
#             Save all exposure SNPs & proxies to be extracted                #
###############################################################################

# Save a list of all exposure SNPs and potential proxies which should be extracted
# as rsid_all.txt; the SNP-outcome associations will be extracted next

all_snps_to_extract <- unique(all_snps)
write(all_snps_to_extract, file.path(out_dir, "rsid_all.txt"))
write.csv(all_proxies, file.path(out_dir, "proxies.csv"))

q('no')