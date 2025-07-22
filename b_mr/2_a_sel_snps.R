###############################################################################
#                                                                             #
#               Select exposure SNPs and proxies locally                      #
#                                                                             #
###############################################################################

###############################################################################
#                                   Set up                                    #
###############################################################################

# Clear the work environment
rm(list=ls())

# Required libraries
## install.packages("remotes")
## remotes::install_github("MRCIEU/TwoSampleMR") 
## remotes::install_github("richardslab/MRutils")
## remotes::install_github("MRCIEU/ieugwasr")
## remotes::install_github("explodecomputer/plinkbinr")
x <- c("dplyr", "tibble", "TwoSampleMR", "MRutils", "ieugwasr", "plinkbinr")
lapply(x, require, character.only = TRUE)

# Set directories
home_dir <- paste0(Sys.getenv("AGE_AT_MENARCHE_DIR"), "working/")
gwas_dir <- paste0(home_dir, "data/b_mr/gwas/")
out_dir <- paste0(home_dir, "/data/b_mr/exposure_dat/")
proxies_dir <- paste0(home_dir, "/data/b_mr/proxies/")
european_reference_panel <- paste0(home_dir, "/data/reference_panels/EUR")
setwd(home_dir)

# Load functions
source("scripts/functions/format_vcf.R")

# Read in exposure GWAS information
opengwas_list <- read.csv(paste0(gwas_dir, "/opengwas_list.csv"))

###############################################################################
#                             Select exposure                                 #
###############################################################################

array_index <- as.numeric((Sys.getenv("SLURM_ARRAY_TASK_ID")))
print(array_index)

phenotype_id <- opengwas_list$phenotype_id[array_index]
print(phenotype_id)

opengwas_id <- opengwas_list$opengwas_id[array_index]
print(opengwas_id)

###############################################################################
#                     Clump instruments and find proxies                      #
###############################################################################

file_path <- paste0(out_dir, phenotype_id)
print(paste0("Phenotype: ", phenotype_id, ", ID: ", opengwas_id))

### Read in and format VCF
print("Formatting VCF")
tmp <- format_vcf(id = opengwas_id, vcf_file_path = gwas_dir,
        exposure_name = phenotype_id)
head(tmp)

### Filter to genome-wide significant SNPs only
tmp_genomewidesig <- filter(tmp, pval.exposure < 5e-08)
rm(tmp) # to reduce memory burden
if(nrow(tmp_genomewidesig)>0){
    genomewidesig_path <- paste0(file_path, "_pval_5e_08_raw.csv")
    write.csv(tmp_genomewidesig, genomewidesig_path, row.names = F)
    print(paste("Genome-wide significant SNPs saved:", genomewidesig_path))
    head(tmp_genomewidesig)
}
    
### Clump
if(nrow(tmp_genomewidesig)>0){
    print("Clumping SNPs")
    # Rename for ld_clump function
    tmp_clumped <- rename(tmp_genomewidesig, pval = pval.exposure, rsid = SNP)
    
    # To run ld_clump through the ieugwas API use the below line however the server can be busy
    # so this command won't always run
    # tmp_clumped <-  ld_clump(tmp_clumped, clump_kb=10000, clump_r2=0.001, pop = "EUR")
    # To run ld_clump locally, use the below line.
    # See: https://mrcieu.github.io/ieugwasr/articles/local_ld.html
    # To download reference panel:
    # wget -nc http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz data/reference_panels/
    tmp_clumped <-  ld_clump(tmp_clumped, clump_kb = 10000, clump_r2 = 0.001,
                      plink_bin = get_plink_exe(), # needs plinkbinr package loaded
                      bfile = european_reference_panel)
    
    # Rename back to TwoSampleMR format
    tmp_clumped$id <- NULL
    tmp_clumped <- rename(tmp_clumped, pval.exposure = pval, SNP = rsid)
    
    clumped_path <- paste0(file_path, "_pval_5e_08_clumped.csv") 
    write.csv(tmp_clumped, clumped_path, row.names = F)
    print(paste("Clumped SNPs saved:", clumped_path))
    head(tmp_clumped)
}

### Find proxies for all exposure SNPs, r2 > 0.8
if(nrow(tmp_genomewidesig)>0){
    proxies <- get_proxies(
    rsids = tmp_clumped$SNP,
    token = Sys.getenv("LDLINK_TOKEN"),
        # request a personal LDlink token here:
        # https://ldlink.nih.gov/?tab=apiaccess
    population = "CEU",
    results_dir = proxies_dir,
    skip_api = FALSE,
    r2_threshold = 0.8
    )
    proxies_path <- paste0(file_path, "_proxies.csv")
    write.csv(proxies, proxies_path, row.names = F)
    print(paste("Possible proxy SNPs saved:", proxies_path))
    head(proxies)
}

q('no')