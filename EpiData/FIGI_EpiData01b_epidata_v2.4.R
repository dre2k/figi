#=============================================================================#
# EPI DATA V2.4
# 03/09/2020
# 
# Incorporate E variables
# Process duplicates, problem samples
# Subset to gxe set
#
#
# ------ NOTES ------ #
# 
# see documentation and update this script
#
#=============================================================================#
library(data.table) 
library(tidyverse)
setwd("~/data/FIGI_samplefile_epi-200309/")
#-----------------------------------------------------------------------------#
# edit colo23 csv file pooledcompassid (only run once)
#-----------------------------------------------------------------------------#

# colo <- read.csv('epi/109colo230.csv', header=T, fileEncoding='latin1', stringsAsFactors= F, na.strings=c('',NA)) %>% 
#   mutate(pooledcompassid = gsub("&", "", pooledcompassid))
# write.csv(colo, file = "epi/109colo230.csv", quote = F, row.names = F)


#-----------------------------------------------------------------------------#
# incorporate epi data
#-----------------------------------------------------------------------------#
rm(list = ls())

# start with samplefile, filter by drops/gxe
# load("~/data/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_190729.RData")
load("~/data/FIGI_EpiData_rdata/FIGI_SampleFile_Summary_200309.RData")

df <- figi_samplefile %>% 
  filter(drop == 0)
any(duplicated(df$pooledcompassid))


# when reading epi data, don't include variables names found in samplefile (except for pooledcompassid)
list_of_variables <- read.table("~/data/FIGI_samplefile_epi-191205/list_of_variables.txt", stringsAsFactors = F) %>% 
  filter(!V1 %in% names(figi_samplefile)[!names(figi_samplefile) %in% 'pooledcompassid']) %>% 
  .[,1]

wrap <- function(x) {
  read.csv(x, colClasses = 'character', fileEncoding='latin1', stringsAsFactors= F, na.strings=c('',NA)) %>% 
    dplyr::select_if(names(.) %in% list_of_variables)
}
tmp <- lapply(list.files("~/data/FIGI_samplefile_epi-200309/epi", full.names = T), wrap)


# find intersection of column names for all epi datasets
# run wrap function again with this set of variable names
list_of_variables_epi <- Reduce(intersect, lapply(tmp, names))

wrap <- function(x) {
  read.csv(x, colClasses = 'character', fileEncoding='latin1', stringsAsFactors= F, na.strings=c('',NA)) %>% 
    dplyr::select_if(names(.) %in% list_of_variables_epi)
}
tmp <- do.call(rbind, lapply(list.files("~/data/FIGI_samplefile_epi-200309/epi", full.names = T), wrap))


# final figi set
figi <- inner_join(df, tmp, by = 'pooledcompassid')
  

any(duplicated(figi$vcfid))
any(duplicated(figi$pooledcompassid))



# Epi data, no genotype information
mismatches <- anti_join(tmp, df, 'pooledcompassid')
table(mismatches$gecco_study)

mismatches <- anti_join(df, tmp, 'pooledcompassid')
table(mismatches$study)

rm(list=setdiff(ls(), c("figi_samplefile", "figi")))
save.image("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_200309.RData")



#-----------------------------------------------------------------------------#
# Sample counts for Yi
#-----------------------------------------------------------------------------#
rm(list = ls())
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_200309.RData")

# gwas set
counts <- data.frame(table(figi$study, figi$outc)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  rename(study = Var1)
write.csv(counts, file = "~/Dropbox/HRC_v2.3_GWAS_Set_counts.csv", quote = F, row.names = F)

# gwas set
gxe <- filter(figi, gxe == 1)
counts <- data.frame(table(gxe$study, gxe$outc)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  rename(study = Var1)
write.csv(counts, file = "~/Dropbox/HRC_v2.3_GxE_Set_counts.csv", quote = F, row.names = F)





#-----------------------------------------------------------------------------#
# Some error checking, accounting for mismatches
# 
#-----------------------------------------------------------------------------#
library(tidyverse)
library(data.table)
library(figifs)
rm(list = ls())

# compare july and december files.. 
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")
july <- figi
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")
december <- figi

# non-union for gwas set
# a single sample changed (CCFR_4) - OK!
non_matches <- anti_join(july, december, by = 'vcfid')
table(non_matches$drop)
table(non_matches$study_gxe)


# gxe set
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190729.RData")
july <- filter(figi, gxe == 1)
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_191205.RData")
december <- filter(figi, gxe == 1)

non_matches <- anti_join(july, december, by = 'vcfid')



#--------------------------------------------#
# quick check asp_ref and aspirin for Flora
#--------------------------------------------#

asp_ref <- figi %>% 
  filter(drop == 0, gxe == 1, EUR_subset == 1, 
         !study_gxe %in% c("ColoCare_1", "ESTHER_VERDI", "GALEON", "MOFFITT", "NFCCR_1"),
         !is.na(asp_ref))
table(asp_ref$asp_ref)

drops <- data.frame(table(asp_ref$study_gxe, asp_ref$outc)) %>% 
  filter(Freq == 0)

asp_ref <- asp_ref %>% 
  filter(!study_gxe %in% drops$Var1)

table(asp_ref$outc)
table(asp_ref$study_gxe, asp_ref$asp_ref, useNA = 'ifany')
table(asp_ref$outc)

table(asp_ref$outc, is.na(asp_ref$aspirin))




asp_ref <- figi %>% 
  filter(drop == 0, gxe == 1, EUR_subset == 1, 
         !study_gxe %in% c("ColoCare_1", "ESTHER_VERDI", "GALEON", "MOFFITT", "NFCCR_1"),
         !is.na(aspirin))

drops <- data.frame(table(asp_ref$study_gxe, asp_ref$outc)) %>% 
  filter(Freq == 0)

asp_ref <- asp_ref %>% 
  filter(!study_gxe %in% drops$Var1)

table(asp_ref$outc)
table(asp_ref$study_gxe, asp_ref$asp_ref, useNA = 'ifany')
table(asp_ref$outc)

table(asp_ref$outc, is.na(asp_ref$aspirin))



