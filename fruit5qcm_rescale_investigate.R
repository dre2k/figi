# explore rescaling for fruit5qcm variable
#
# per Jim, I'll use average sex/study specific IQR as a scaling factor, then run meta-analysis
# this is done overall (which makes sense)
#
# Q why is this better than mean 0 sd 1? because those are better suited for inherently continuous variables, not ones based on quartile medians
# plus the distribution might be skewed, not normal


# CAN LIKELY DELETE AFTER YOU FIGURE IT OUT

# -------------- input data ------------------ #
exposure = 'fruit5qcm'
hrc_version = 'v3.0'
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = '/media/work/gwis_test/fruit5qcm/'

esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset) %>% 
  mutate(fruit5qcm = as.numeric(fruit5qcm))



# ----- calculate IQR by study/sex ------ #
out <- input_data %>% 
  group_by(study_gxe, sex) %>% 
  mutate(fruit5qcm.iqr = IQR(fruit5qcm)) %>% 
  ungroup()

hist(out$fruit5qcm.iqr)

fruit5qcm.iqr.median = median(out$fruit5qcm.iqr)

out <- out %>% 
  mutate(fruit5qcm_scaled = fruit5qcm / fruit5qcm.iqr.median)

summary(out$fruit5qcm.iqr)

# test dachs_1
dachs1 <- filter(out, study_gxe == "DACHS_1")
table(dachs1$fruit5qcm_scaled)


# output forest plot
tmp_forest(data_epi = out, 
           exposure = 'fruit5qcm_scaled', covariates = covariates, hrc_version = hrc_version, path = "~/Dropbox/Working/", categorical = F)








# --------------- 1) overall ------------ #

out <- input_data %>% 
  mutate(fruit5qcm = scale(fruit5qcm))

dachs1 <- filter(out, study_gxe == "DACHS_1")

table(dachs1$fruit5qcm)



# --------------- 2) by study  ------------ #

out <- input_data %>% 
  group_by(study_gxe) %>% 
  mutate(fruit5qcm = scale(fruit5qcm)) %>% ungroup()

dachs1 <- filter(out, study_gxe == "DACHS_1")

table(dachs1$fruit5qcm)


# output forest plot
tmp_forest(data_epi = out, 
           exposure = 'fruit5qcm', covariates = covariates, hrc_version = hrc_version, path = "~/Dropbox/Working/", categorical = F)



