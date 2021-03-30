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

posthoc_sig_wrapper <- function(filename, output_dir) {
  
  # check if data.frame is empty (no rows)
  tmp <- readRDS(glue(output_dir, filename)) %>% 
    mutate(method = gsub(".rds", "", filename))
  
  if (grepl("manhattan", filename)) {
    out <- tmp %>% 
      dplyr::rename(Pval = P) %>% 
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval, method) %>% 
      dplyr::arrange(Chromosome, Location)
  } else {
    out <- tmp %>% 
      dplyr::rename(Pval = step2p) %>% 
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval, method) %>% 
      dplyr::arrange(Chromosome, Location)
  }
}


#----------------------------------------------------------------------#
# Note: 
# - GxE 1df suggestive hits (5e-6), which can be a lot of markers
#   Use LD clumped @ 5e-6 to identify top SNPs
#----------------------------------------------------------------------#

# Filter GxE, 2df, 3df suggestive results to top SNPs only
gxe_clump <- fread(glue("/media/work/gwis/clump_combined/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGxE_ldclump.clumped"), data.table = F) %>% pull(SNP)
gxe_snps <- readRDS(glue("/media/work/gwis/posthoc/{exposure}/manhattan_chiSqGxE_{exposure}_df.rds")) %>%
  filter(SNP %in% gxe_clump)
saveRDS(gxe_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/manhattan_chiSqGxE_", exposure, "_clump_df.rds"))

control_clump <- fread(glue("/media/work/gwis/clump_combined/FIGI_{hrc_version}_gxeset_{exposure}_chiSqControl_ldclump.clumped"), data.table = F) %>% pull(SNP)
control_snps <- readRDS(glue("/media/work/gwis/posthoc/{exposure}/manhattan_chiSqControl_{exposure}_df.rds")) %>%
  filter(SNP %in% control_clump)
saveRDS(control_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/manhattan_chiSqControl_", exposure, "_clump_df.rds"))

case_clump <- fread(glue("/media/work/gwis/clump_combined/FIGI_{hrc_version}_gxeset_{exposure}_chiSqCase_ldclump.clumped"), data.table = F) %>% pull(SNP)
case_snps <- readRDS(glue("/media/work/gwis/posthoc/{exposure}/manhattan_chiSqCase_{exposure}_df.rds")) %>%
  filter(SNP %in% case_clump)
saveRDS(case_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/manhattan_chiSqCase_", exposure, "_clump_df.rds"))

ge_clump <- fread(glue("/media/work/gwis/clump_combined/FIGI_{hrc_version}_gxeset_{exposure}_chiSqGE_ldclump.clumped"), data.table = F) %>% pull(SNP)
ge_snps <- readRDS(glue("/media/work/gwis/posthoc/{exposure}/manhattan_chiSqGE_{exposure}_df.rds")) %>%
  filter(SNP %in% ge_clump)
saveRDS(ge_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/manhattan_chiSqGE_", exposure, "_clump_df.rds"))



twodf_clump <- fread(glue("/media/work/gwis/clump_combined/FIGI_{hrc_version}_gxeset_{exposure}_chiSq2df_no_gwas_ldclump.clumped"), data.table = F) %>% pull(SNP)
twodf_snps <- readRDS(glue("/media/work/gwis/posthoc/{exposure}/manhattan_chiSq2df_{exposure}_no_gwas_df.rds")) %>%
  filter(SNP %in% twodf_clump)
saveRDS(twodf_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/manhattan_chiSq2df_", exposure, "_no_gwas_clump_df.rds"))

threedf_clump <- fread(glue("/media/work/gwis/clump_combined/FIGI_{hrc_version}_gxeset_{exposure}_chiSq3df_no_gwas_ldclump.clumped"), data.table = F) %>% pull(SNP)
threedf_snps <- readRDS(paste0("/media/work/gwis/posthoc/", exposure, "/manhattan_chiSq3df_", exposure, "_no_gwas_df.rds")) %>%
  filter(SNP %in% threedf_clump)
saveRDS(threedf_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/manhattan_chiSq3df_", exposure, "_no_gwas_clump_df.rds"))


filelist <- read.table(text = glue("
manhattan_chiSqGxE_{exposure}_clump_df.rds
manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds
manhattan_chiSq3df_{exposure}_no_gwas_no_ge_clump_df.rds
manhattan_chiSqControl_{exposure}_clump_df.rds
manhattan_chiSqGE_{exposure}_df.rds
twostep_wht_chiSqEDGE_{exposure}_df.rds 
twostep_wht_chiSqGE_{exposure}_gwas_step1_no_gwas_df.rds
twostep_wht_chiSqEDGE_{exposure}_expectation_hybrid_df.rds
twostep_wht_chiSqGE_{exposure}_no_gwas_df.rds
twostep_wht_chiSqEDGE_{exposure}_gwas_step1_df.rds
twostep_wht_chiSqG_{exposure}_df.rds
twostep_wht_chiSqEDGE_{exposure}_gwas_step1_no_gwas_df.rds
twostep_wht_chiSqG_{exposure}_expectation_hybrid_df.rds
twostep_wht_chiSqEDGE_{exposure}_no_gwas_df.rds
twostep_wht_chiSqG_{exposure}_gwas_step1_df.rds
twostep_wht_chiSqGE_{exposure}_df.rds
twostep_wht_chiSqG_{exposure}_gwas_step1_no_gwas_df.rds
twostep_wht_chiSqGE_{exposure}_expectation_hybrid_df.rds
twostep_wht_chiSqG_{exposure}_no_gwas_df.rds
twostep_wht_chiSqGE_{exposure}_gwas_step1_df.rds
manhattan_chiSq2df_{exposure}_functional_subset_control_no_gwas_df.rds
manhattan_chiSqGxE_{exposure}_functional_subset_control_no_gwas_df.rds
manhattan_chiSq2df_{exposure}_functional_subset_pooled_no_gwas_df.rds
manhattan_chiSqGxE_{exposure}_functional_subset_pooled_no_gwas_df.rds
manhattan_chiSq2df_{exposure}_functional_subset_tumor_no_gwas_df.rds
manhattan_chiSqGxE_{exposure}_functional_subset_tumor_no_gwas_df.rds
manhattan_chiSq3df_{exposure}_functional_subset_control_no_gwas_no_ge_df.rds
twostep_wht_chiSqEDGE_{exposure}_functional_subset_control_no_gwas_df.rds
manhattan_chiSq3df_{exposure}_functional_subset_pooled_no_gwas_no_ge_df.rds
twostep_wht_chiSqEDGE_{exposure}_functional_subset_pooled_no_gwas_df.rds
manhattan_chiSq3df_{exposure}_functional_subset_tumor_no_gwas_no_ge_df.rds 
twostep_wht_chiSqEDGE_{exposure}_functional_subset_tumor_no_gwas_df.rds
manhattan_chiSqCase_{exposure}_functional_subset_control_no_gwas_df.rds
twostep_wht_chiSqGE_{exposure}_functional_subset_control_no_gwas_df.rds
manhattan_chiSqCase_{exposure}_functional_subset_pooled_no_gwas_df.rds
twostep_wht_chiSqGE_{exposure}_functional_subset_pooled_no_gwas_df.rds
manhattan_chiSqCase_{exposure}_functional_subset_tumor_no_gwas_df.rds
twostep_wht_chiSqGE_{exposure}_functional_subset_tumor_no_gwas_df.rds
manhattan_chiSqG_{exposure}_functional_subset_control_no_gwas_df.rds
twostep_wht_chiSqG_{exposure}_functional_subset_control_no_gwas_df.rds
manhattan_chiSqG_{exposure}_functional_subset_pooled_no_gwas_df.rds
twostep_wht_chiSqG_{exposure}_functional_subset_pooled_no_gwas_df.rds
manhattan_chiSqG_{exposure}_functional_subset_tumor_no_gwas_df.rds
twostep_wht_chiSqG_{exposure}_functional_subset_tumor_no_gwas_df.rds
twostep_wht_chiSqEDGE_{exposure}_expectation_hybrid_no_gwas_df.rds
twostep_wht_chiSqG_{exposure}_expectation_hybrid_no_gwas_df.rds
twostep_wht_chiSqGE_{exposure}_expectation_hybrid_no_gwas_df.rds"), header = F)



filelist <- read.table(text = glue("
manhattan_chiSqGxE_{exposure}_clump_df.rds
manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds
manhattan_chiSq3df_{exposure}_no_gwas_clump_df.rds
manhattan_chiSqCase_{exposure}_clump_df.rds
manhattan_chiSqControl_{exposure}_clump_df.rds"), header = F)







# calcium_totqc2
filelist <- read.table(text = glue("
manhattan_chiSqGxE_{exposure}_clump_df.rds
manhattan_chiSqGE_calcium_totqc2_clump_df.rds
manhattan_chiSqControl_calcium_totqc2_clump_df.rds
manhattan_chiSqCase_calcium_totqc2_clump_df.rds
manhattan_chiSq2df_{exposure}_no_gwas_clump_df.rds
manhattan_chiSq3df_{exposure}_no_gwas_clump_df.rds
twostep_wht_chiSqGE_calcium_totqc2_expectation_hybrid_df.rds"), header = F)



# ------ this might be necessary depending on the exposure ------- #
# 2df clumped if needed (clump because often the peaks have many SNPs in LD)
# clump, then remove gwas hits to decreased number of SNPs
# exclude_gwas <- fread("~/data/Annotations/gwas_140_chr_bp_ref_alt_plink.tags", header = F)$V1
# exclude_gwas <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt", header = F)$V1
# twodf_clump <- fread(paste0("/media/work/gwis/clump_combined/FIGI_", hrc_version, "_gxeset_", exposure, "_chiSq2df_ldclump.clumped")) %>%
#   filter(!SNP %in% exclude_gwas) %>%
#   pull(SNP)
# twodf_snps <- readRDS(paste0("/media/work/gwis/posthoc/", exposure, "/significant_results_dataframe_chiSq2df_", exposure, ".rds")) %>%
#   filter(SNP %in% twodf_clump) %>%
#   rename(Pval_2df = Pval) # this is just to play nice with older scripts where i DIDN"T clump chiSq2df stats
# saveRDS(twodf_snps, file = paste0("/media/work/gwis/posthoc/", exposure, "/significant_results_dataframe_chiSq2df_", exposure, "_clump.rds"))

# --------------------------------------------------------#
# read the significant result data.frame, output file for binarydosage getsnp function

# manhattan results
tmp1 <- as.character(filelist$V1)
tmp2 <- lapply(tmp1, function(x) posthoc_sig_wrapper(x, output_dir = output_dir) )
out <- do.call(rbind, tmp2)

# add Rsq and alt allele frequency information for these hits (for convenience when creating posthoc results report)
info <- readRDS("/media/work/FIGI_RsqEstimate_chrALL.rds") %>% 
  dplyr::filter(id %in% out$SNP) %>% 
  dplyr::select(-maf, -SNP)
out2 <- inner_join(out, info, by = c("SNP" = "id"))

saveRDS(out2, file = glue(output_dir, "gwis_sig_results_input_{exposure}.rds"), version = 2)







# ------ alcoholc_moderate ------- #
out_extra <- c("9:97251034:C:T", 9, 97251034, "C", "T", 64937, 27733 , 0.05, "gwas2")

out2 <- rbind(out, out_extra) %>% 
  filter(method == "gwas2")
saveRDS(out2, file = glue(output_dir, "gwis_sig_results_additional_input_{exposure}.rds"), version = 2)


# ------ calcium_totqc2 ------- #
# exclude E|G among controls (no need to extract those)
# also exclude corresponding Case-only results for that particular locus

ca1 <- filter(out2, method == "manhattan_chiSqCase_calcium_totqc2_df" & Location == 114550461)

ca2 <- out2 %>% 
  dplyr::filter(method != "manhattan_chiSqGE_calcium_totqc2_df", 
                method != "manhattan_chiSqCase_calcium_totqc2_df")

ca_out <- rbind(ca2, ca1)

saveRDS(ca_out, file = glue(output_dir, "gwis_sig_results_input_{exposure}.rds"), version = 2)



# ------ calcium_dietqc2 ------- #
# exclude E|G among controls (no need to extract those)
# also exclude corresponding Case-only results for that particular locus

ca1 <- out2 %>% 
  dplyr::filter(method != "manhattan_chiSqCase_calcium_dietqc2_df",
                method != "manhattan_chiSqGE_calcium_dietqc2_df", 
                method != "manhattan_chiSqControl_calcium_dietqc2_df")

saveRDS(ca1, file = glue(output_dir, "gwis_sig_results_input_{exposure}.rds"), version = 2)



# ------ folate_totqc2 ------- #
# MTHFR locus, might want to consider those in LD as well... 
# rs1801133 (1:11856378:G:A)

out_extra <- c("1:11856378:G:A", 1, 11856378, "G", "A", 58464, 27655, 0.05, "additional_locus")

out2 <- rbind(out, out_extra) %>% 
  filter(method == "additional_locus")

saveRDS(out2, file = glue(output_dir, "gwis_sig_results_additional_input_{exposure}.rds"), version = 2)




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
