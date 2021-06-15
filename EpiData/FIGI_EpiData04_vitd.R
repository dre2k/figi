#=============================================================================#
# explore vitamin D distribution etc 
#=============================================================================#
library(tidyverse)
library(data.table)
library(stargazer)
library(lmtest)
library(broom)
library(figifs)
library(kableExtra)
library(qqman)
rm(list = ls())

# FIGI gxe set (gxe == 1, EUR_subset == 1)
figi_old <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_analysis_data_glm.rds")
figi <- readRDS("/media/work/gwis/data/FIGI_EpiData/FIGI_v2.3_gxeset_analysis_data_glm.rds")

load("/home/rak/data/FIGI_vitaminD_epi-200527/pooled_vitD_25OH.RData")
load("/home/rak/data/FIGI_vitaminD_epi-200527/vitd_crp_UKB.Rdata")


table(vitD$blood_month, vitD$study, useNA = 'ifany')

# can i get country information from epic studies
epic <- filter(figi, grepl("EPIC", study_gxe ))
table(epic$study_site, useNA = 'ifany')


hpfs <- filter(figi, grepl("HPFS", study_gxe ))
#-----------------------------------------------------------------------------#
# atbc, epic, hpfs, mccs, mec, nhs, whi
#-----------------------------------------------------------------------------#
# counts are correct, missing N=21 due to EUR_subset != 1
figiD <- inner_join(figi, vitD, 'pooledcompassid') %>% 
  mutate(group = "others")


figiD_outlier <- filter(figiD, vitd25 < 1000)

qqnorm(figiD_outlier$vitd25)
qqline(figiD_outlier$vitd25)

qqPlot(figiD_outlier$vitd25)

qqPlot(log(figiD_outlier$vitd25))




# 
# 
# 
# 

spain <- c("Asturias", "Granada", "Murcia", "Navarra", "San Sebastian")
netherlands <- c("Bilthoven", "Utrecht")
uk <- c("Cambridge", "UK General Population", "UK Health Conscious")
italy <- c("Florence", "Naples", "Turin", "Varese")
greece <- c("Greece")
germany <- c("Heidelberg", "Potsdam")
france <- c("Northeast France", "Northwest France", "South coast France", "South France")
croatia <- c("Ragusa")
switzerland <- c("Umea")


epic <- filter(figiD, grepl("EPIC", study_gxe )) %>% 
  mutate(study_country = case_when(study_site %in% spain ~ "Spain", 
                                   study_site %in% netherlands ~ "Netherlands", 
                                   study_site %in% uk ~ "UK", 
                                   study_site %in% italy ~ "Italy", 
                                   study_site %in% greece ~ "Greece", 
                                   study_site %in% germany ~ "Germany", 
                                   study_site %in% france ~ "France", 
                                   study_site %in% croatia ~ "Croatia", 
                                   study_site %in% switzerland ~ "Switzerland"))
table(epic$study_site, useNA = 'ifany')
table(epic$study_site, epic$study_country, useNA = 'ifany')

table(epic$study_country, epic$outc)

#-----------------------------------------------------------------------------#
# ukb 
# check if 25-hydroxy (is it consistent with other dataset) 
#-----------------------------------------------------------------------------#
bio.uk <- mutate(bio.uk, f.eid = as.character(f.eid))
ukb <- filter(figi, study_gxe == "UKB_1")
ukbD <- inner_join(ukb, bio.uk, by = c("compassid" = "f.eid")) %>% 
  mutate(vitd25 = f.30890.0.0, 
         crp = f.30710.0.0, 
         group = 'ukb')

hist(bio.uk$f.30890.0.0)
hist(ukbD$vitd25)

#-----------------------------------------------------------------------------#
# 
#-----------------------------------------------------------------------------#
out <- bind_rows(figiD[, c('pooledcompassid', 'group', 'vitd25', 'studyname')], 
                 ukbD[, c('pooledcompassid', 'group', 'vitd25', 'studyname')]) %>% 
  mutate(group = as.factor(group)) 

out_mean <- out %>% 
  group_by(group) %>%
  dplyr::summarize(group.mean = mean(vitd25, na.rm = T))

out_mean <- out %>% 
  group_by(studyname) %>%
  dplyr::summarize(group.mean = mean(vitd25, na.rm = T))


summary(out$vitd25)

ggplot(out, aes(group, vitd25)) + 
  geom_boxplot()

ggplot(out, aes(studyname, vitd25)) +
  geom_boxplot()

# two outliers - remove if vitd25 > 1000
# (both from EPIC study)

out <- filter(out, vitd25 < 1000)
ggplot(out, aes(group, vitd25)) + 
  geom_boxplot()

ggplot(out, aes(studyname, vitd25)) + 
  geom_boxplot()

ggplot(out, aes(vitd25, color = group, fill = group)) + 
  geom_histogram(alpha = 0.5, position = 'identity') + 
  geom_vline(data=out_mean, aes(xintercept=group.mean, color=group), linetype="dashed")

ggplot(out, aes(vitd25, color = studyname, fill = studyname)) + 
  geom_histogram(alpha = 0.5, position = 'identity') + 
  geom_vline(data=out_mean, aes(xintercept=group.mean, color=studyname), linetype="dashed")










#-----------------------------------------------------------------------------#
# create dataset for gxescanr... . . . 
#-----------------------------------------------------------------------------#

# (just for reference since this is how i created all the other datasets)
wrap <- function(d, exposure, is_e_categorical = T, min_cell_size = 0, vars_to_exclude = c("energytot_imp"), vars_to_include = c(), studies_to_exclude) {
  cov <- format_data_glm(d, exposure, is_e_categorical, min_cell_size, vars_to_exclude, vars_to_include, eur_only = T) %>% 
    filter(!study_gxe %in% studies_to_exclude) %>% 
    mutate(study_gxe = fct_drop(study_gxe))
  saveRDS(cov, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds"), version = 2)
  
  cov_gxescan <- format_data_gxescan(cov, exposure)
  saveRDS(cov_gxescan, file = paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_gxescan.rds"), version = 2)
}



# merge nonukb studies
# call variable vitd25
# (counts are correct, missing N=21 due to EUR_subset != 1)
figiD <- inner_join(figi, vitD, 'pooledcompassid') %>% 
  mutate(group = "others")
figiD_outlier <- filter(figiD, vitd25 < 1000)


# ukbiobank
bio.uk <- mutate(bio.uk, f.eid = as.character(f.eid))
ukb <- filter(figi, study_gxe == "UKB_1")
ukbD <- inner_join(ukb, bio.uk, by = c("compassid" = "f.eid")) %>% 
  mutate(vitd25 = f.30890.0.0, 
         crp = f.30710.0.0, 
         group = 'ukb')

variable_list <- c( "vcfid", "outcome", "age_ref_imp", "sex", "study_gxe", "pc1", "pc2", "pc3", "vitd25")

out <- bind_rows(figiD_outlier[, variable_list], 
                 ukbD[, variable_list]) %>% 
  mutate(study_gxe = fct_drop(study_gxe))
saveRDS(out, file = "/media/work/gwis/results/input/FIGI_v2.3_gxeset_vitd25_basic_covars_glm.rds")

out_gxescan <- format_data_gxescan(out, 'vitd25')
saveRDS(out_gxescan, file = "/media/work/gwis/results/input/FIGI_v2.3_gxeset_vitd25_basic_covars_gxescan.rds", version = 2)





# let's look at batch effects for WHI, EPIC, and BioUK (aliquot??) and distribution of vitD




tmp <- vitD %>%
  filter(study == "WHI") %>% 
  group_by(whi_testver) %>% 
  summarise(vitd_mean = mean(as.numeric(vitd25)))


tmp <- vitD %>%
  filter(study == "EPIC") %>% 
  group_by(epic_batch_vitd) %>% 
  summarise(vitd_mean = mean(as.numeric(vitd25)))


tmp <- vitD %>%
  filter(study == "WHI") %>% 
  mutate(whi_testver = factor(whi_testver, labels = c('2', '3','14')))


ggplot(tmp, aes(whi_testver, vitd25, group = whi_testver)) + 
  geom_boxplot()

tmp <- vitD %>%
  filter(study == "EPIC", 
         vitd25 < 1000) %>% 
  mutate(epic_batch_vitd = factor(epic_batch_vitd))


ggplot(tmp, aes(epic_batch_vitd, vitd25, group = epic_batch_vitd)) + 
  geom_boxplot()


tmp <- bio.uk %>%
  filter(f.eid %in% figi$compassid) %>% 
  mutate(f.30892.0.0 = factor(f.30892.0.0))

ggplot(tmp, aes(f.30892.0.0, f.30890.0.0, group = f.30892.0.0)) + 
  geom_boxplot()




#-----------------------------------------------------------------------------#
# association testing
#-----------------------------------------------------------------------------#

# out_glm <- bind_rows(figiD, ukbD) %>% 
#   filter(vitd25 < 1000) %>% 
#   mutate(vitd25_log = log(vitd25), # natural log
#          vitd25_sqrt = sqrt(vitd25)) # sqrt
# 
# qqPlot(out_glm$vitd25)
# qqPlot(out_glm$vitd25_log)
# qqPlot(out_glm$vitd25_sqrt)
# 
# 
# summary(glm(outcome ~ vitd25 + sex + age_ref_imp + study_gxe, data = out_glm, family = 'binomial'))
# summary(glm(outcome ~ vitd25_log + sex + age_ref_imp + study_gxe, data = out_glm, family = 'binomial'))
# summary(glm(outcome ~ vitd25_sqrt + sex + age_ref_imp + study_gxe, data = out_glm, family = 'binomial'))

## -------- create like a summary table of means, sds etc for each study. then also the same for histograms and barplots, and qq plots. Then we can get a good idea 
## 
## 
## 
## 
## 
## 


# ------------------ vitamin D main effects ----------------------- #
# have to recreate the analysis while creating updated variables.. 

figi <- readRDS("/media/work/gwis/data/FIGI_EpiData/FIGI_v2.3_gxeset_analysis_data_glm.rds")
load("/home/rak/data/FIGI_vitaminD_epi-200527/pooled_vitD_25OH.RData")
load("/home/rak/data/FIGI_vitaminD_epi-200527/vitd_crp_UKB.Rdata")


# create dataset
tmp1 <- inner_join(figi, dplyr::select(vitD, -study), 'pooledcompassid') %>% 
  mutate(group = "others")
tmp1_outlier <- filter(tmp1, vitd25 < 1000)

# ukbiobank
tmp2 <- mutate(bio.uk, f.eid = as.character(f.eid))
tmp3 <- filter(figi, study_gxe == "UKB_1")
tmp4 <- inner_join(tmp3, tmp2, by = c("compassid" = "f.eid")) %>% 
  mutate(vitd25 = f.30890.0.0, 
         crp = f.30710.0.0, 
         group = 'ukb')

# out file
# 
# EPIC countries
spain <- c("Asturias", "Granada", "Murcia", "Navarra", "San Sebastian")
netherlands <- c("Bilthoven", "Utrecht")
uk <- c("Cambridge", "UK General Population", "UK Health Conscious")
italy <- c("Florence", "Naples", "Turin", "Varese", "Ragusa")
greece <- c("Greece")
germany <- c("Heidelberg", "Potsdam")
france <- c("Northeast France", "Northwest France", "South coast France", "South France")
switzerland <- c("Umea")


out <- bind_rows(tmp1_outlier, tmp4) %>% 
  mutate(study_country = case_when(study_site %in% spain ~ "Spain", 
                                   study_site %in% netherlands ~ "Netherlands", 
                                   study_site %in% uk ~ "UK", 
                                   study_site %in% italy ~ "Italy", 
                                   study_site %in% greece ~ "Greece", 
                                   study_site %in% germany ~ "Germany", 
                                   study_site %in% france ~ "France", 
                                   study_site %in% switzerland ~ "Switzerland", 
                                   study %in% c("NHS", "HPFS", "WHI") ~ "USA", 
                                   study == "UKB" ~ "UK"), 
         study_gxe_tmp1 = case_when(study == 'EPIC' ~ study_country,
                                       TRUE ~ study_gxe), 
         study_gxe_tmp2 = case_when(study_gxe == "WHI_1" & whi_testver == 2 ~ "WHI_1_2", 
                                    study_gxe == "WHI_1" & whi_testver == 3 ~ "WHI_1_3",
                                    study_gxe == "WHI_1" & whi_testver == 14 ~ "WHI_1_14",
                                    study_gxe == "WHI_2" & whi_testver == 2 ~ "WHI_2_2", 
                                    study_gxe == "WHI_2" & whi_testver == 3 ~ "WHI_2_3",
                                    study_gxe == "WHI_2" & whi_testver == 14 ~ "WHI_2_14",
                                    study_gxe == "WHI_3" & whi_testver == 2 ~ "WHI_3_2", 
                                    study_gxe == "WHI_3" & whi_testver == 3 ~ "WHI_3_3",
                                    study_gxe == "WHI_3" & whi_testver == 14 ~ "WHI_3_14", 
                                    TRUE ~ study_gxe_tmp1))

# ugh
table(out$study_gxe, out$whi_testver)
table(out$study_gxe_tmp2, out$outc)

model1 <- glm(outcome ~ vitd25 + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
summary(model1)


model2 <- glm(outcome ~ vitd25 + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe_tmp1, data = out, family = 'binomial')
summary(model2)


model3 <- glm(outcome ~ vitd25 + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe_tmp2, data = out, family = 'binomial')
summary(model3)



