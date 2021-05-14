#=============================================================================#
# instead of filtering GWAS snps by LD, just remove SNP BP +- n_flank regions
#
# same as above but using GIANT consortium meta-analysis results
# might be best to remove entire regions to be honest... 
#=============================================================================#


# 1) clump the meta analysis results to get only top SNPs in each peak. then you can remove based on flanking region
clumped <- fread("/media/work/gwis_test/data/Meta-analysis_Wood_et_al_UKBiobank_2018_clump.clumped") %>% 
  pull(SNP)

clumped <- fread("/scratch2/andreeki/gwis_test/data/Meta-analysis_Wood_et_al_UKBiobank_2018_clump.clumped") %>% 
  pull(SNP)


# 2) get regions for all 7000 items (machine might crash..)
gxe <- qread("FIGI_v2.3_gxeset_height10_cohort_basic_covars_gxescan_results.qs")

test <- function(data, snp, n_flank = 500000) {
  snp_sep = unlist(strsplit(snp, split = ":"))
  chr = as.numeric(snp_sep[1])
  bp = as.numeric(snp_sep[2])

  out <- data %>%  
    filter(Chromosome == chr & between(Location, bp-n_flank, bp+n_flank)) %>% 
    pull(SNP2)
  
  return(out)
}


# Map over list of gwas hits (chr:bp)
tmp <- test(gxe, clumped[1])
tmp <- unlist(map(clumped[1:2], ~ test(gxe, .x)))








# Testing things out
gwas <- fread("~/Dropbox/FIGI/FIGI_code/code/Annotations/update_210226/conditioning_snps_v20200930.tsv")
giant <- fread("~/Meta-analysis_Wood_et_al_UKBiobank_2018.txt")
gxe <- qread("/media/work/gwis/results/aspirin/processed/FIGI_v2.4_gxeset_aspirin_basic_covars_gxescan_results.qs")


# write function, then map over vector of gwas chr:bp
gwas_list <- gwas %>% 
  mutate(SNP2 = paste0(CHR, ":", POS)) %>% 
  pull(SNP2) %>% unique(.)

test <- function(data, snp, n_flank = 500000) {
  
  snp_sep = unlist(strsplit(snp, split = ":"))
  chr = as.numeric(snp_sep[1])
  bp = as.numeric(snp_sep[2])
  
  out <- data %>%  
    filter(Chromosome == chr & between(Location, bp-n_flank, bp+n_flank)) %>% 
    pull(SNP2)
  
  return(out)
}

# tmp <- gxe %>% 
#   filter(Chromosome == 1 & between(Location, 22503282-500000, 22503282+500000))

# Map over list of gwas hits (chr:bp)
tmp <- test(gxe, "1:22503282")
tmp <- unlist(map(gwas_list[1:2], ~ test(gxe, .x)))



# --------------------------------------------------------------------------- #
# execute, save as RDS
# use this vector to filter SNPs from processed gxe results
# --------------------------------------------------------------------------- #

tmp <- unique(unlist(map(clumped, ~ test(gxe, .x))))
saveRDS(tmp, file = "~/Dropbox/FIGI/FIGI_code/code/Annotations/update_210226/conditioning_snps_v20200930_filter_GWAS_SNPS.rds")
