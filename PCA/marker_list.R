# =========================================================================== #
# PCA revisit
# =========================================================================== #

# marker list for PCA calculation

x <- fread("~/figi_gxe_2500samples_v2_nold_maf1_prune.prune.in", stringsAsFactors = F, header = F)


# this is a list of SNPs with rsq > 0.8
rsq_filter <- readRDS("~/data/Rsq_Estimate/FIGI_RsqEstimate_chrALL.rds")


rsq <- rsq_filter %>% 
  mutate(snp = sub('^([^:]+:[^:]+).*', '\\1', id)) %>% 
  filter(snp %in% x$V1)

# these are the SNPs MAF > 1% and also Rsq > 0.8. randomly sample 30,000 of them
# (also removed high LD regions)
set.seed(2021)
out <- sample(rsq$snp, 30000)

write.table(out, file = "~/figi_markerlist_update.txt", quote = F, row.names = F, col.names = F)
