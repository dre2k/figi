test <- readRDS("/media/work/gwis_test/data/FIGI_v2.3_aaf_gxeset_eur_noUKB_chr22.rds")


test2 <- readRDS("/media/work/gwis_test/data/snps_with_allele_freq_greater_than_5_noukb.rds")


genome <- 1:22
genome
chromos <- 1:22
chr= 22
bdose <- readRDS(glue("/media/work/gwis_test/bdose_info/FIGI_snpid_fix_chr{chr}.rds"))
aaf <- readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_aaf_gxeset_eur_noUKB_chr{chr}.rds"))
out <- data.frame(snp = bdose$SNPs,
                  aaf = aaf)
head(out)
out <- data.frame(snp = bdose$SNPs,
                  aaf_noukb = aaf)
out <- data.frame(snp = bdose$SNPs,
                  aaf_noukb = aaf) %>%
  filter(aaf_noukb >= 0.05)
out <- data.frame(snp = bdose$SNPs,
                  aaf_noukb = aaf) %>%
  filter(aaf_noukb >= 0.05) %>%
  pull(snp)
out <- data.frame(snp = bdose$SNPs,
                  aaf_noukb = aaf) %>%
  filter(aaf_noukb >= 0.05)
head(out)
out <- data.frame(bdose$SNPs,
                  aaf_noukb = aaf) %>%
  filter(aaf_noukb >= 0.05)
head(out)
out <- data.frame(bdose$SNPs,
                  aaf_noukb = aaf) %>%
  filter(aaf_noukb >= 0.05) %>%
  pull(SNPID)

wrap <- function(chr) {
  bdose <- readRDS(glue("/media/work/gwis_test/bdose_info/FIGI_snpid_fix_chr{chr}.rds"))
  aaf <- readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_aaf_gxeset_eur_noUKB_chr{chr}.rds"))
  out <- data.frame(bdose$SNPs,
                    aaf_noukb = aaf) %>%
    filter(aaf_noukb >= 0.01) %>%
    pull(SNPID)
  return(out)
}
aaf_out <- do.call(c, map(chromos, ~ wrap(.x)))
head(aaf_out)
saveRDS(aaf_out, "/media/work/gwis_test/data/snps_with_allele_freq_greater_than_0.01_noukb.rds")


