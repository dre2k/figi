dd <- fread("~/toAndre2.csv")




#Note that all of the discrepant records are indicated as a case on the phys Act data,
# and as a control in the older gwas data.  Last column of attached is the indicator of high/low physical activity exposure based on the 8.75 met hours cutoff.







# v2.3
figi_23 <- readRDS("/media/work/gwis_test/data/FIGI_v2.3_GWAS.rds")
figi_31 <- readRDS("/media/work/gwis_test/data/FIGI_v3.1_GWAS.rds") 

table(figi_23$study_gxe)

tmp <- filter(figi_23, vcfid %in% dd$vcfid)

table(tmp$study_gxe)
table(tmp$outcome)
tmp$FID

tmp <- filter(figi_31, vcfid %in% dd$vcfid)
table(tmp$outcome)


tmp <- filter(figi_31, grepl("1630_HIPFX_550D_4", vcfid))
tmp <- filter(figi_31, grepl("6150376017_R02C01", ))
tmp <- filter(figi_31, grepl("6150376017_R02C01", vcfid))
tmp <- filter(figi_31, grepl("6150376017_R02C01", vcfid))
tmp <- filter(figi_31, grepl("6150376017_R02C01", vcfid))
tmp <- filter(figi_31, grepl("6150376017_R02C01", vcfid))


tmp <- filter(figi_31, netcdfid %in% dd$V1)

table(figi_23$study_gxe)

table(figi_23$outcome, figi_23$V2)
