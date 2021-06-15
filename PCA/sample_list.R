# =========================================================================== #
# PCA revisit
# =========================================================================== #

# SNPs for PCA should be pruned, and large LD regions should be excluded
# 
# see:
# Price, Alkes L., Michael E. Weale, Nick Patterson, Simon R. Myers, Anna C. Need, Kevin V. Shianna, Dongliang Ge, et al. 2008. “Long-Range LD Can Confound Genome Scans in Admixed Populations.” American Journal of Human Genetics.
# 
# Weale, Michael E. 2010. “Quality Control for Genome-Wide Association Studies.” Methods in Molecular Biology  628: 341–72.
# 
# Wellcome Trust Case Control Consortium. 2007. “Genome-Wide Association Study of 14,000 Cases of Seven Common Diseases and 3,000 Shared Controls.” Nature 447 (7145): 661–78.
# 
# 
# (although i honestly don't recall anyone saying anything even after I presented on how I performed PCA)
# in any case - let's explore whether this will make any difference in the results, starting witn NSAIDs



# ------- Steps: ------- #
# - sample 2500 individuals from FIGI GxE set
# - obtained all SNPs, prune, and remove LD regions (see https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD))
# - make sure those SNPs are high Rsq (maybe MAF > 1% as well)
# - sample 30,000 from these, then extract from VCF files for all GxE set samples (takes a while)
# - THEN perform PCA, obtain eigenvectors, perform QC, and fit models with known hits for aspirin/asp_ref to check


# --------------------------------------------------------------------------- #
# begin 
# --------------------------------------------------------------------------- #

figi <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_v3.0_gxeset_analysis_data_glm.rds"))

set.seed(2021)
sample_list <- sample(figi$vcfid, size = 2500, replace = F)
out <- data.frame(id1 = sample_list, id2 = sample_list)

write.table(out, file = "~/figi_vcfid_sample_5000.txt", quote = F, row.names = F)

check <- filter(figi, vcfid %in% sample_list)
table(check$study_gxe)




# --------------------------------------------------------------------------- #
# output vcfids for gxe set (to extract the 30,000 SNPs for PCA)
# --------------------------------------------------------------------------- #

figi <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_v3.0_gxeset_analysis_data_glm.rds"))
out <- data.frame(id1 = figi$vcfid, id2 = figi$vcfid)
write.table(out, file = "~/figi_vcfid_gxeset.txt", quote = F, row.names = F)

gwas <- readRDS("~/data/FIGI_EpiData_rdata/")

gwas <- readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_GWAS.rds")) %>% 
  mutate(id1 = vcfid, id2 = vcfid) %>% 
  dplyr::select(id1, id2)
write.table(gwas, file = "~/figi_vcfid_gwas.txt", quote = F, row.names = F, sep = '\t')


