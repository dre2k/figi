#=============================================================================#
# GxEScanR Results pre-process
# Updated 05/27/2020
# Updated 08/18/2020
# 
# Inputs
# - gxescanr output
# - Rsq estimates (recalculated using John's code)
# 
# Tasks
# - filter SNPs Rsq > 0.8
# - calculate EDGE statistics 
# - calculate 3df statistics
#
# Output
# - rds results file
# - p values for LD clumping by chromosome ('obsolete')
# - p values for LD clumping, all chromosomes combined
# - p values for LocusZoom plots
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
library(glue)
library(qs)
rm(list = ls())

# parse arguments (if running from command line)
# args <- commandArgs(trailingOnly=T)
# exposure <- args[1] # ex: asp_ref
# hrc_version <- args[2] # ex: v2.4
exposure = 'aspirin'
hrc_version = 'v2.4'
wdir <- "/media/work/gwis/results/"
output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")

#-----------------------------------------------------------------------------#
# Read results, filter by Rsq
#-----------------------------------------------------------------------------#
# Rsq estimate based on alt allele probability, maf > 0.001, rsq >= 0.8
# (you can find this function in my fork of John's BinaryDosage R package)
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")
qsave(rsq_filter, file = "~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.qs")
rsq_filter <- qread("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.qs")

gxe_all <- do.call(rbind, lapply(list.files(path = paste0(wdir, exposure), full.names = T, pattern = paste0("FIGI_", hrc_version, "_gxeset_", exposure)), fread, stringsAsFactors = F, data.table = F)) %>%
  mutate(SNP2 = paste0(Chromosome, ":", Location), 
         chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE, 
         chiSqEDGE_p = pchisq(chiSqEDGE, df = 2, lower.tail = F),
         chiSqG_p = pchisq(chiSqG, df = 1, lower.tail = F), 
         chiSqGxE_p = pchisq(chiSqGxE, df = 1, lower.tail = F), 
         chiSq2df_p = pchisq(chiSq2df, df = 2, lower.tail = F), 
         chiSq3df_p = pchisq(chiSq3df, df = 3, lower.tail = F), 
         chiSqGE_p = pchisq(chiSqGE, df = 1, lower.tail = F), 
         chiSqCase_p = pchisq(chiSqCase, df = 1, lower.tail = F), 
         chiSqControl_p = pchisq(chiSqControl, df = 1, lower.tail = F))

# if necessary .. 
# gxe_all <- mutate(gxe_all, SNP = paste0(Chromosome, ":", Location, ":", Reference, ":", Alternate))
gxe <- gxe_all %>%
  filter(SNP %in% rsq_filter$id)

# create directory, save rds file
system(paste0("mkdir -p ", paste0(wdir, exposure, "/processed/")))

qsave(gxe, file = paste0(wdir, exposure, paste0("/processed/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_gxescan_results.qs")))
# saveRDS(gxe, file = paste0(wdir, exposure, paste0("/processed/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_gxescan_results.rds")), version = 2)



#-----------------------------------------------------------------------------#
# Read-in processed file, for convenience when re-running
#-----------------------------------------------------------------------------#
gxe <- qread(paste0(wdir, exposure, paste0("/processed/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_gxescan_results.qs")))

#-----------------------------------------------------------------------------#
# Comments:
# 
# - for two-step methods, LD clumping based on GxE statistics is biased, omit
#   clumping on D|G seems more 'fair'. Ultimately, we decided against genome-
#   wide clumping. Instead, will use expectation-based hybrid method
# 
# - when creating manhattan plots of 2df/3df results, it's useful to remove
#   significant D|G and E|G hits to check if there are loci driven more by GxE
#   Use LD clumping to identify tags + SNPs in LD for D|G and E|G, and remove
#   from plots
# - in retrospect, I'd rather visually inspect all peaks. Limit removal of 
#   peaks to the 140 GWAS loci from Jeroen paper
#   
# - clump markers genome-wide to easily identify the top hit in each region. 
#   this makes outputting results table easier
#   
# - the folder '/media/work/gwis' contains all results + posthoc analysis
#-----------------------------------------------------------------------------#

### ------- output p values by chromosome ------- ###
# Don't do this anymore
# for(chr in 1:22) {
#   out <- gxe %>%
#     # dplyr::mutate(P = pchisq(chiSqG, df = 1, lower.tail = F)) %>%
#     dplyr::rename(P = paste0(statistic, "_p")) %>%
#     dplyr::filter(Chromosome == chr) %>%
#     dplyr::select(SNP, P)
#   write.table(out, file = paste0("/media/work/gwis/clump_bychr/", exposure, "/", "FIGI_", hrc_version, "_gxeset_", exposure,  "_basic_covars_gxescan_chr", chr, "_chiSqG_ldclump.txt"), quote = F, row.names = F, sep = '\t')
# }

### ------- output p values, single file ------- ###
ldclump_wrapper <- function(statistic) {
  out <- gxe %>%
    dplyr::rename(P = paste0(statistic, "_p")) %>%
    dplyr::filter(P < 5e-4) %>% 
    dplyr::select(SNP, P)
  
  write.table(out, file = paste0("/media/work/gwis/clump_combined/FIGI_",
                                 hrc_version, "_gxeset_", exposure,
                                 paste0("_", statistic, "_ldclump.txt")),
              quote = F, row.names = F, sep = '\t')
}

ldclump_wrapper('chiSqG')
ldclump_wrapper('chiSqGxE')
ldclump_wrapper('chiSqControl')
ldclump_wrapper('chiSqCase')
ldclump_wrapper('chiSqGE')

# also want to clump 2df/3df but after removing gwas loci (to identify new hits)
# need to create 'exclude_gwas' object first.. 
ldclump_wrapper_nogwas <- function(statistic, df) {
  out <- gxe %>%
    dplyr::rename(P = paste0(statistic, "_p")) %>%
    dplyr::filter(P < 5e-4) %>% 
    dplyr::filter(!SNP2 %in% exclude_gwas) %>% 
    dplyr::select(SNP, P)
  
  write.table(out, file = paste0("/media/work/gwis/clump_combined/FIGI_",
                                 hrc_version, "_gxeset_", exposure,
                                 paste0("_", statistic, "_no_gwas_ldclump.txt")),
              quote = F, row.names = F, sep = '\t')
}

exclude_gwas <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt") %>%
  mutate(SNP2 = paste0(Chr, ":", Pos)) %>%
  pull(SNP2)


exclude_gwas <- fread("~/data/Annotations/gwas_200_ld_annotation_feb2021.txt") %>%
  mutate(SNP2 = paste0(Chr, ":", Pos)) %>%
  pull(SNP2)


ldclump_wrapper_nogwas('chiSq2df', df = 2)
ldclump_wrapper_nogwas('chiSq3df', df = 3)

#-----------------------------------------------------------------------------#
# Output p-values for locuszoom plots
# 
# check: 
# - how does locuszoom handle duplicate positions?? 
#-----------------------------------------------------------------------------#
locuszoom_wrapper <- function(statistic, df) {
  locuszoom <- gxe %>%
    dplyr::rename(`P-value` = paste0(statistic, "_p")) %>%
    dplyr::mutate(MarkerName = paste0("chr", Chromosome, ":", Location)) %>%
    dplyr::select(MarkerName, `P-value`)

  write.table(locuszoom, file = paste0("/media/work/gwis/locuszoom/FIGI_",
                                       hrc_version, "_gxeset_", exposure,
                                       paste0("_basic_covars_gxescan_", statistic, "_locuszoom.txt")),
              quote = F, row.names = F, sep = "\t")
}

locuszoom_wrapper('chiSqG', df = 1)
locuszoom_wrapper('chiSqGxE', df = 1)
locuszoom_wrapper('chiSqControl', df = 1)
locuszoom_wrapper('chiSqCase', df = 1)
locuszoom_wrapper('chiSq2df', df = 2)
locuszoom_wrapper('chiSq3df', df = 3)



#-----------------------------------------------------------------------------#
# two-step expectation hybrid ---- 
#-----------------------------------------------------------------------------#
# ------ output list of SNPs from each bin (expectation based) to extract dosages
# exposure subset for step 1 statistics

output_dir_eh <- glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/')
system(paste0("mkdir -p ", output_dir_eh))

twostep_eh_snps(gxe, 'chiSqG', output_dir = output_dir_eh)
twostep_eh_snps(gxe, 'chiSqGE', output_dir = output_dir_eh)
twostep_eh_snps(gxe, 'chiSqEDGE', output_dir = output_dir_eh)

# gecco meta analysis for main effects step 1 statistics
gwas_results <- readRDS("/media/work/gwis/results/gecco/MarginalMeta_gecco_HRC_EUR_only.rds")
gxe_twostep_gwas_step1 <- gxe %>%
  dplyr::select(-betaG, -chiSqG, -chiSqEDGE, -chiSq3df) %>%
  inner_join(gwas_results, 'SNP') %>%
  mutate(chiSqEDGE = chiSqG + chiSqGE,
         chiSq3df = chiSqG + chiSqGxE + chiSqGE)

twostep_eh_snps(gxe_twostep_gwas_step1, 'chiSqG', step1_source = "gecco", output_dir = output_dir_eh)
