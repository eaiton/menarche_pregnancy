############################################################################
#                  Function to add proxies to outcome data                 #
############################################################################

# Locally generates outcome data including proxy SNPS, if available

# requires library(tidyr)

# outcome_dat is the full outcome dataset, including all proxy SNPs
# proxies is the list of proxies retrieved using MRutils::get_proxies in
# 2_a_sel_snps.R

generate_outdat_with_proxies <- function(exposure_dat, outcome_dat, outcome_name,
                                         proxies){
  
  ### Set up
  i_outcome <- outcome_name
  # Filter specific outcome dataset to available exposure SNPs
  tmp_outcome_dat <- outcome_dat %>%
    filter(SNP %in% exposure_dat$SNP, outcome == i_outcome)
  
  ### Remove any duplicated SNPs
  tmp_outcome_dat <- tmp_outcome_dat[!duplicated(tmp_outcome_dat$SNP), ]
  
  ### Check whether outcome is available
  if(dim(tmp_outcome_dat)[1] == 0){
    print(paste("Outcome", outcome_name, "not available"))
    
    # Return data with 0 rows if outcome was not available
    tmp_outcome_dat
    
    ### Check whether proxies needed
  } else if (dim(tmp_outcome_dat)[1] == dim(exposure_dat)[1]){
    print(paste("No proxies needed for", outcome_name))
    
    # Return full outcome dataset if no proxies needed
    tmp_outcome_dat
    
  } else {
    ### Identify proxies
    
    # Specify all available SNPs
    outcome_snps <- expand.grid(
      SNP = c(exposure_dat$SNP), outcome = unique(outcome_dat$outcome))
    # Left join outcome_dat, NAs if missing
    outcome_snps_available <-
      left_join(outcome_snps, outcome_dat, by = c("SNP", "outcome")) %>%
      # Add indicator column for NAs
      mutate(missing = ifelse(is.na(effect_allele.outcome),
                              TRUE, FALSE))
    # Deduplicated list of all SNPs which we need a proxy for
    need_proxies <- outcome_snps_available %>%
      filter(missing == TRUE) %>%
      distinct(SNP)
    # Count of proxies needed by outcome
    table(outcome_snps_available$outcome[outcome_snps_available$missing == TRUE])
    # List of all proxies
    proxy_snps <- unique(proxies$rsid)
    # Proxies without query SNPs
    proxies_merge_diff <- subset(proxies, query_rsid!=rsid)
    proxies_merge_diff$SNP<-proxies_merge_diff$rsid
    
    outcome_dat_proxies_tmp <-subset(outcome_dat,
                                     outcome_dat$outcome== i_outcome)
    
    # All SNPs we need to proxy for this outcome
    outcome_snps_tmp <- subset(outcome_snps, outcome_snps$outcome==i_outcome)
    # Available SNPs in outcome data:
    outcome_snps_available <-
      left_join(outcome_snps_tmp, outcome_dat, by = c("SNP", "outcome")) %>%
      # Add indicator column for NAs
      mutate(missing = ifelse(is.na(effect_allele.outcome),
                              TRUE, FALSE))
    # SNPs not available for our outcome in outcome data, for which we need proxies:
    need_proxies <- outcome_snps_available %>%
      filter(missing == TRUE) %>%
      distinct(SNP)
    
    if(dim(need_proxies)[1] > 0){
      
      print(paste0("Searching for proxies for ", dim(need_proxies)[1], " SNPs."))
      
      # Append proxies information -
      # For every SNP needing a proxy ('SNP'),
      # create row with information about the proxy ('rsid' & 'SNP.y' columns)
      tmp1 <- merge(x=need_proxies, y=proxies_merge_diff, by.x="SNP", by.y="query_rsid", all.x=TRUE)
      
      # For each proxy SNP ('rsid'), append the extracted proxy-outcome GWAS data
      tmp2 <-merge(x=tmp1, y=outcome_dat_proxies_tmp, by.x="rsid", by.y="SNP", all.y=TRUE)
      # this leaves NAs where any proxy SNPs could not be extracted
      # since they were not available in outcome data
      
      # Exclude palindromic proxies:
      tmp2 <-subset(tmp2, !(effect_allele.outcome == "A" & other_allele.outcome == "T" |
                              effect_allele.outcome == "T" & other_allele.outcome == "A") )
      tmp2 <-subset(tmp2, !(effect_allele.outcome == "C" & other_allele.outcome == "G" |
                              effect_allele.outcome == "G" & other_allele.outcome == "C") )
      
      # Select the proxy SNP ('rsid') in highest LD with the query SNP ('SNP') needing a proxy -
      set.seed(67898) # to ensure slice_sample() is reproducible
      tmp3 <- tmp2 %>%
        group_by(SNP) %>%
        # Remove any proxies with missing outcome data
        drop_na() %>%
        # Select max R2 value, keep ties if several have same R2
        slice_max(R2, with_ties = TRUE) %>%
        # Select random proxy SNP if several have same R2
        slice_sample(n = 1)
      
    }
    
    # Stop searching for proxies if no SNP-outcome data available
    if( (dim(need_proxies)[1] > 0) & (dim(tmp3)[1] == 0) ){
      
      print("No proxy-outcome associations available for any SNPs.")
      
    } else {
      
      # Add column specifying proxy effect allele
      # checking proxy outcome data and using original SNP alleles
      # using query SNP alleles but checking they're the right way round using proxy outcome data
      #'Correlated alleles' format:
      # query effect allele = proxy effect allele, query other allele = proxy other allele
      tmp4 <- tmp3
      tmp4$effect_allele.proxy <- apply(tmp4, 1, function(row) {
        
        # If outcome effect allele same as proxy effect allele, treat as if using query effect allele -
        if (row["effect_allele.outcome"] == substr(row["Correlated_Alleles"], 3, 3)) {
          return(substr(row["Correlated_Alleles"], 1, 1))
          # if outcome effect allele same as proxy other allele, treat as if using query other allele -
        } else if (row["effect_allele.outcome"] == substr(row["Correlated_Alleles"], 7, 7)) {
          return(substr(row["Correlated_Alleles"], 5, 5))
        } else {
          return("default_value")
        }
      })
      # Same again but for other allele -
      tmp4$other_allele.proxy <- apply(tmp4, 1, function(row) {
        if (row["other_allele.outcome"] == substr(row["Correlated_Alleles"], 3, 3)) {
          return(substr(row["Correlated_Alleles"], 1, 1))
        } else if (row["other_allele.outcome"] == substr(row["Correlated_Alleles"], 7, 7)) {
          return(substr(row["Correlated_Alleles"], 5, 5))
        } else {
          return("default_value")
        }
      })
      
      # Treat proxy alleles as the outcome alleles
      tmp5 <- tmp4
      tmp5$effect_allele.outcome<-tmp5$effect_allele.proxy
      tmp5$other_allele.outcome<-tmp5$other_allele.proxy
      tmp5$data_source.outcome <- NA
      
      # Create final dataframe which has replaced missing SNPs with their proxies
      outcome_dat_cols <- colnames(outcome_dat)
      tmp5 <- tmp5 %>% select(all_of(outcome_dat_cols))
      
      # Return message with number of SNPs which could be proxied
      print(paste0(dim(tmp5)[1], " proxies identified."))
      
      tmp_outcome_dat <- rbind(tmp_outcome_dat, tmp5)
    }
    
  }
  
  # Return full outcome dataset with added proxies
  tmp_outcome_dat
  
}