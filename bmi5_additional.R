bmi5 <- qread("/media/work/gwis_test/bmi5/results/processed/FIGI_v2.3_gxeset_bmi5_basic_covars_gxescan_results.qs")

snplist <- c("11:120531547:C:G",
             "11:46574459:A:T",
             "11:47328213:C:T",
             "15:33122966:C:T",
             "2:10961902:C:T",
             "2:217410559:A:G",
             "3:100037644:A:G",
             "6:138766021:C:T",
             "9:136601110:G:C",
             "10:114274269:T:C",
             "12:50177324:T:A",
             "15:32987718:G:A",
             "15:33122966:C:T",
             "20:6584196:C:G",
             "5:134486618:G:A",
             "15:33122966:C:T")

snplist <- c("11:120531547",
             "11:46574459",
             "11:47328213",
             "15:33122966",
             "2:10961902",
             "2:217410559",
             "3:100037644",
             "6:138766021",
             "9:136601110",
             "10:114274269",
             "12:50177324",
             "15:32987718",
             "15:33122966",
             "20:6584196",
             "5:134486618",
             "15:33122966")

unique(snplist)

bmi5_short <- dplyr::select(bmi5, SNP, Chromosome, Location, Reference, Alternate,  Subjects, Cases, contains("_p")) %>% 
  mutate(SNP2 = paste0(Chromosome, ":", Location)) %>% 
  filter(SNP2 %in% snplist) 
write.csv(bmi5_short, file = "~/Dropbox/Working/figi_bmi_pvalues.csv")
rm(bmi5)  





# females

bmi5 <- qread("/media/work/gwis_test/bmi5_female/results/processed/FIGI_v2.3_gxeset_bmi5_female_basic_covars_gxescan_results.qs")

snplist <- c("1:243472801",
             "11:18966592",
             "11:47328213",
             "19:51341455",
             "19:56597188",
             "2:20972616",
             "2:21095479",
             "2:36580097",
             "2:45771123",
             "4:84609954",
             "5:91530803",
             "6:119181125",
             "7:5144105")


bmi5_short <- dplyr::select(bmi5, SNP, Chromosome, Location, Reference, Alternate,  Subjects, Cases, contains("_p")) %>% 
  mutate(SNP2 = paste0(Chromosome, ":", Location)) %>% 
  filter(SNP2 %in% snplist) 

write.csv(bmi5_short, file = "~/Dropbox/Working/figi_bmi_female_pvalues.csv")








# males

bmi5 <- qread("/media/work/gwis_test/bmi5_male/results/processed/FIGI_v2.3_gxeset_bmi5_male_basic_covars_gxescan_results.qs")

snplist <- c("19:12483407",
             "22:49189233",
             "3:115329531",
             "4:185473157",
             "6:160627010",
             "8:133357436",
             "8:15169150")


bmi5_short <- dplyr::select(bmi5, SNP, Chromosome, Location, Reference, Alternate,  Subjects, Cases, contains("_p")) %>% 
  mutate(SNP2 = paste0(Chromosome, ":", Location)) %>% 
  filter(SNP2 %in% snplist) 

write.csv(bmi5_short, file = "~/Dropbox/Working/figi_bmi_male_pvalues.csv")







# 
# 
# 
# hrc <- fread("/media/work/ReferencePanels/HRC.r1-1.GRCh37.wgs.mac5.sites.tab")
# 
# hrc_short <- hrc %>% 
#   dplyr::select(`#CHROM`, POS, ID, REF, ALT) %>% 
#   dplyr::mutate(SNP = paste(`#CHROM`, POS, REF, ALT, sep = ':'))
# rm(hrc)
# 
# 
# head(bmi5_short)
# 
# out  <- inner_join(bmi5_short, hrc_short[, c("SNP", "ID")], 'SNP')
# 
# # tmp <- head(out) %>% 
# #   mutate(chiSqG_p = pchisq(chiSqG, df = 1, lower.tail = F), 
# #          seG = betaG / sqrt(chiSqG))
# # head(tmp)
# 
# # have his group deal with the NaNs. 
# qc <- out %>% 
#   filter(!is.finite(seGxE))
# 
# 
# write.csv(out, file = "~/tmp.csv", quote = F, row.names = F)
# saveRDS(out, file = "~/tmp.rds")
