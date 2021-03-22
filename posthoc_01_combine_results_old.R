#=============================================================================#
# post-hoc analysis of FIGI GWIS findings
# 03/18/2020
#
# Notes:
# - take output from rmarkdown report 
#   (data.frames RDS files of significant and suggestive results)
#   and rbind them 
# - P values from 1/2/3df tests, P values from step 2 statistic for 2step method
#
# - Columns:
#   "SNP"
#   "Chromosome"
#   "Location"
#   "Reference"
#   "Alternate"
#   "Subjects"
#   "Cases"
#   "Pval"
#   "method"
#
# Create 'input' files - serve as input for getSNPvalues (HPC)
# Input files also useful to get list of SNPs for further posthoc analysis
#
# treat expectation based findings separately (in another script perhaps)
#=============================================================================#
library(tidyverse)
library(data.table)
library(glue)

posthoc_sig_wrapper <- function(filename) {

  # check if data.frame is empty (no rows)
  tmp <- readRDS(paste0("/media/work/gwis/posthoc/", exposure, "/", filename))
  
  if(is.null(dim(tmp)) | dim(tmp)[1] == 0) {
    out <- data.frame()
  } else {
    if(grepl("twostep_wht", filename)) {
      # if twostep finding, output step2p pvalue (GxE 1DF)
      out <- tmp %>% 
        dplyr::rename(Pval = step2p) %>% 
        dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval) %>% 
        mutate(method = gsub("significant_results_dataframe_|.rds", "", filename))
    } else if(grepl("2df", filename)) {
      # if 2df finding, output 2df statistic
      out <- tmp %>%
        dplyr::rename(Pval = Pval_2df) %>% 
        dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval) %>% 
        mutate(method = gsub("significant_results_dataframe_|.rds", "", filename))
    } else if(grepl("3df", filename)) {
      # if 3df finding, output 3df statistic
      out <- tmp %>%
        dplyr::rename(Pval = Pval_3df) %>% 
        dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval) %>% 
        mutate(method = gsub("significant_results_dataframe_|.rds", "", filename))
    } else {
      # otherwise - statistic (GxE 1df, case only - 'Pval')
      out <- tmp %>% 
        dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval) %>% 
        mutate(method = gsub("significant_results_dataframe_|.rds", "", filename))
    }
  } 
}


#----------------------------------------------------------------------#
# 
#
# Note: 
# - GxE 1df suggestive hits (5e-6), which can be a lot of markers
#   Use LD clumped @ 5e-6 to identify top SNPs
#----------------------------------------------------------------------#

# GxE 1df clumped findings
gxe_clump <- fread(paste0("/media/work/gwis/clump_combined/FIGI_", hrc_version, "_gxeset_", exposure, "_chiSqGxE_ldclump.clumped"))$SNP
gxe_snps <- readRDS(paste0("/media/work/gwis/posthoc/", exposure, "/significant_results_dataframe_chiSqGxE_", exposure, ".rds")) %>%
  filter(SNP %in% gxe_clump)
saveRDS(gxe_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/significant_results_dataframe_chiSqGxE_", exposure, "_clump.rds"))

# 2df clumped if needed (clump because often the peaks have many SNPs in LD)
# clump, then remove gwas hits to decreased number of SNPs
exclude_gwas <- fread("~/data/Annotations/gwas_140_chr_bp_ref_alt_plink.tags", header = F)$V1
exclude_gwas <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt", header = F)$V1
twodf_clump <- fread(paste0("/media/work/gwis/clump_combined/FIGI_", hrc_version, "_gxeset_", exposure, "_chiSq2df_ldclump.clumped")) %>%
  filter(!SNP %in% exclude_gwas) %>%
  pull(SNP)
twodf_snps <- readRDS(paste0("/media/work/gwis/posthoc/", exposure, "/significant_results_dataframe_chiSq2df_", exposure, ".rds")) %>%
  filter(SNP %in% twodf_clump) %>%
  rename(Pval_2df = Pval) # this is just to play nice with older scripts where i DIDN"T clump chiSq2df stats
saveRDS(twodf_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/significant_results_dataframe_chiSq2df_", exposure, "_clump.rds"))


# --------------------------------------------------------#
# edit this part as needed 

#  TO DO RIGHT NOW (4/16/2020) - ADD OTHER FILES YOUD LIKE TO GET DOSAGES FOR
filelist <- read.table(text = glue("
significant_results_dataframe_chiSq2df_{exposure}_clump.rds
significant_results_dataframe_chiSq3df_{exposure}_no_gwas_no_ge.rds
significant_results_dataframe_chiSqCase_{exposure}.rds
significant_results_dataframe_chiSqGxE_{exposure}_clump.rds
significant_results_dataframe_twostep_wht_chiSqEDGE_{exposure}.rds
significant_results_dataframe_twostep_wht_chiSqEDGE_{exposure}_chiSqG_ld_clump.rds
significant_results_dataframe_twostep_wht_chiSqG_{exposure}.rds
significant_results_dataframe_twostep_wht_chiSqG_{exposure}_chiSqG_ld_clump.rds
significant_results_dataframe_twostep_wht_chiSqGE_{exposure}.rds
significant_results_dataframe_twostep_wht_chiSqGE_{exposure}_chiSqG_ld_clump.rds"), header = F)

# exclude expectation based results
filelist <- read.table(text = glue("
significant_results_dataframe_chiSqGxE_{exposure}_clump.rds
significant_results_dataframe_chiSq2df_{exposure}_no_gwas_no_marginal.rds
significant_results_dataframe_chiSq3df_{exposure}_no_gwas_no_marginal_no_ge.rds
significant_results_dataframe_chiSqCase_{exposure}.rds"), header = F)



# GxE clumped
# case only
# 2df no gwas, clumped
# 3df no gwas, clumped
# two step wht test for three statistic (mostly empty)
# two step ld clumped
# ADD EXPECTATION BASED SOMETIME
filelist <- read.table(text = glue("
significant_results_dataframe_chiSqGxE_{exposure}_ld_clump_no_gwas.rds
significant_results_dataframe_chiSqCase_{exposure}.rds
significant_results_dataframe_chiSq2df_{exposure}_ld_clump_no_gwas.rds
significant_results_dataframe_chiSq3df_{exposure}_ld_clump_no_gwas.rds
significant_results_dataframe_twostep_wht_chiSqEDGE_{exposure}.rds
significant_results_dataframe_twostep_wht_chiSqEDGE_{exposure}_chiSqG_ld_clump.rds
significant_results_dataframe_twostep_wht_chiSqG_{exposure}.rds
significant_results_dataframe_twostep_wht_chiSqG_{exposure}_chiSqG_ld_clump.rds
significant_results_dataframe_twostep_wht_chiSqGE_{exposure}.rds
significant_results_dataframe_twostep_wht_chiSqGE_{exposure}_chiSqG_ld_clump.rds"), header = F)

x <- readRDS("/media/work/gwis/posthoc/aspirin/significant_results_dataframe_chiSq2df_aspirin_ld_clump_no_gwas.rds")

# --------------------------------------------------------#

tmp1 <- as.character(filelist$V1)
tmp2 <- lapply(tmp1, posthoc_sig_wrapper)
out <- do.call(rbind, tmp2) %>% 
  mutate(Pval = formatC(as.numeric(Pval), format = "e", digits = 2))

# add Rsq and alt allele frequency information for these hits
info <- readRDS("/media/work/gwis/FIGI_RsqEstimate_chrALL.rds") %>% 
  dplyr::filter(id %in% out$SNP) %>% 
  dplyr::select(-maf, -SNP)
out2 <- inner_join(out, info, by = c("SNP" = "id"))

# save as data.frame
# use to extract allelic dosages from binarydosage files
saveRDS(out2, file = paste0("/media/work/gwis/posthoc/gwis_sig_results_input_", exposure, ".rds"), version = 2)

#-----------------------------------------------------------------------------#
# one-off modifications of the input file 
#-----------------------------------------------------------------------------#

# # ---------- diab ---------- #
# diab <- readRDS("/media/work/tmp/posthoc/gwis_sig_results_input_diab.rds")
# diab_edge_expectation <- readRDS("~/Dropbox/FIGI/Results/diab/files/significant_results_twostep_wht_chiSqEDGE_diab_expectation_based.rds") %>% 
#   filter(SNP %in% c("1:50981205:G:A", "8:118185025:G:A", "11:61621194:G:A", "13:47191972:G:A")) %>% 
#   rename(Pval = step2p) %>% 
#   dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval) %>% 
#   mutate(method = gsub("significant_results_|.rds", "", "significant_results_twostep_wht_chiSqEDGE_diab_expectation_based.rds"))
# 
# diab_out <- rbind(diab, diab_edge_expectation)
# saveRDS(diab_out, file = paste0("/media/work/tmp/posthoc/gwis_sig_results_input_", 'diab', ".rds"), version = 2)
