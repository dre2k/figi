bmi5 <- qread("/media/work/gwis_test/bmi5/results/processed/FIGI_v2.3_gxeset_bmi5_basic_covars_gxescan_results.qs")

bmi5_short <- dplyr::select(bmi5, SNP, Chromosome, Location, Reference, Alternate,  Subjects, Cases, betaG, chiSqG, betaGxE, chiSqGxE) 
rm(bmi5)  

hrc <- fread("/media/work/ReferencePanels/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")


hrc_short <- hrc %>% 
  dplyr::select(`#CHROM`, POS, ID, REF, ALT) %>% 
  dplyr::mutate(SNP = paste(`#CHROM`, POS, REF, ALT, sep = ':'))
rm(hrc)


head(bmi5_short)

out  <- inner_join(bmi5_short, hrc_short[, c("SNP", "ID")], 'SNP') %>% 
  mutate(seG = betaG / sqrt(chiSqG), 
         seGxE = betaGxE / sqrt(chiSqGxE)) %>% 
  dplyr::select(-betaG, -chiSqG, -seG)

# tmp <- head(out) %>% 
#   mutate(chiSqG_p = pchisq(chiSqG, df = 1, lower.tail = F), 
#          seG = betaG / sqrt(chiSqG))
# head(tmp)

# have his group deal with the NaNs. 
qc <- out %>% 
  filter(!is.finite(seGxE))


write.csv(out, file = "~/tmp.csv", quote = F, row.names = F)
saveRDS(out, file = "~/tmp.rds")
