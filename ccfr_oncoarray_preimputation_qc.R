#=====================================================#
# 07/14/2021
#
# try to reimpute but also run HWE on the thinging
# (multi ethnic samples)
# 
# per Jeroen email, run on combined case/control EUR only (self-reported probably?)
# p < 1e-4
#=====================================================# 

sample_list <- fread("~/Dropbox/Working/ccfr_oncoarray_reimputation_20210710/ccfr_oncoarray_hrc_20170828_UW.fam")
ccfr <- fread("~/ccfr_oncoarray_qc_working/101ccfr0_201014.csv")
# load("~/data/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_200309.RData")

# find intersection betwen figi_samplefile and sample_list
out <- ccfr %>% 
  filter(compassid %in% sample_list$V1)

head(ccfr$compassid)
table(out$race_self, out$outc, useNA = 'ifany')
any(duplicated(out$compassid))
table(out$study_site)



# non-matching samples 
non_match <- anti_join(sample_list, out, by = c("V1"= "compassid"))
write.table(non_match, file = "~/Dropbox/ccfr_oncoarray_preimputation_qc_101ccfr_non_matches.txt", quote = F, row.names= F, col.names = F)

# output list of self-reported "White", then use to filter HWE and submit to michican imputation server after running wrayner? 
white <- filter(out, race_self == "White")
write.table(white[, c('compassid', 'compassid')], file = "~/Dropbox/ccfr_oncoarray_white_974.txt", quote = F, row.names = F, col.names = F, sep = '\t')


