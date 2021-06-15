#=============================================================================#
# diabetes GxE additional analyses
# adjusted by other known covariates
#
# don't confuse this with e main effects, this is fitting interaction models
# so i think the best bet would be simply to report estimates in stargazer
#=============================================================================#
library(tidyverse)
library(data.table)
library(figifs)
library(stargazer)
rm(list = ls())

env <- 'diab'

vcfid_list <- x <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", env, "_basic_covars_glm.rds"))[,'vcfid']

geno <- readRDS(paste0("/media/work/tmp/posthoc/gwis_sig_results_output_", env, ".rds")) %>% 
  dplyr::select(-contains('p0'), 
                -contains('p1'),
                -contains('p2'))

figi <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_results_exposure_main_effects.rds") %>% 
  rename_all(tolower) %>%
  rename(EUR_subset = eur_subset) %>% 
  inner_join(geno, 'vcfid') %>% 
  filter(vcfid %in% vcfid_list) %>% 
  mutate(redmeatqc2 = factor(redmeatqc2),
         vigmodlns = as.numeric(vigmodlns))

covars_1 <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
covars_2 <- sort(c(covars_1, 'bmi', 'smk_ever', 'redmeatqc2'))
covars_3 <- sort(c(covars_1, 'bmi', 'smk_ever', 'redmeatqc2', 'vigmodlns'))






snp <- "X10.114813041.C.T_dose"
snp <- "1.71040166"



gxe_glm_wrapper <- function(snp) {
  
  # ----- for diab only, running on subset of SNps, and all you have is chromosome and bp.. 
  snp <- str_subset(names(figi), snp)
  
  # create formulas
  formula_interaction_1 <- as.formula(paste0("outcome ~ ", env, " * ", snp, " + ", paste0(covars_1, collapse = " + ")))
  formula_interaction_2 <- as.formula(paste0("outcome ~ ", env, " * ", snp, " + ", paste0(covars_2, collapse = " + ")))
  formula_interaction_3 <- as.formula(paste0("outcome ~ ", env, " * ", snp, " + ", paste0(covars_3, collapse = " + ")))
  
  # run glm 
  x1 <- glm(formula_interaction_1, data = figi, family = 'binomial')
  x2 <- glm(formula_interaction_2, data = figi, family = 'binomial')
  x3 <- glm(formula_interaction_3, data = figi, family = 'binomial')
  
  coefs <- list(exp(coef(x1)),
                exp(coef(x2)),
                exp(coef(x3)))
  
  snp_clean <- gsub("\\.", "\\_", snp)
  snp_clean <- gsub("X", "chr", snp_clean)
  snp_clean <- gsub("\\_dose", "", snp_clean)
  
  # output html output from stargazer
  out <- capture.output(stargazer(x1, x2, x3, title = paste0(snp_clean, " x T2D with additional adjustment covariates"), align = T, type = 'html', ci=TRUE, ci.level=0.95, omit = c("pc", "study_gxe"), keep.stat = "n", column.labels=c("Covar Set 1","Covar Set 2","Covar Set 3"), star.cutoffs = c(0.05, 0.01, 0.001), column.sep.width = '10pt', coef=coefs, p.auto = F, notes = c("PC and Study/Platform estimates omitted from table", paste0("Covar Set 1 - ", gsub("_", "-", paste0(covars_1, collapse = ","))), paste0("Covar Set 2 - ", gsub("_", "-", paste0(covars_2, collapse = ","))), paste0("Covar Set 3 - ", gsub("_", "-", paste0(covars_3, collapse = ","))))))
  
  cat(paste(out, collapse = "\n"), "\n", file = paste0("/media/work/tmp/results_report/gxe_glm_", snp_clean, "_x_", env, "_additional_covars.html"), append = F)
  # cat(paste(out, collapse = "\n"), "\n", file = paste0("~/Dropbox/gxe_glm_", snp_clean, "_x_", env, "_additional_covars.html"), append = F)
  
}



#
#
# do some processing since you're only doing on subset and you have the basepair locations ... 

# str_subset(names(figi), "1.50981205" )

snp_list <- c("1.71040166", 
              "3.185510884",
              "6.20688121",
              "7.28189411",
              "8.118185025",
              "10.114754071",
              "10.114784926",
              "12.4384844",
              "16.53811788",
              "20.6442961")


for(x in snp_list) {
  gxe_glm_wrapper(x)
}


out <- c()
for(x in snp_list) {
  h <- str_subset(names(figi), x)
  out <- c(out, h)
}






file_list <- c(unlist(list.files("/media/work/tmp/results_report/", pattern = '^gxe_glm_chr.*diab_additional_covars.html?', full.names = T)))





#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
### for Diab... a little more complicated. you need GxE results for likelihood ratio test for 2df/3df
### for 2df... ok. 
### for 3DF .. hm.. not sure. I might have to do it the way that jim mentions it.. 
### - it does work, the estimates are SLIGHTLY different because of the way John implemented the 
### G|E models, confirm with him

# let's test it out right now

# 2df on 13:47191972
str_subset(names(figi), "13.47191972")
env<-'diab'

glm_func <- glm(outcome ~ X13.47191972.G.A_dose * diab + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
glm_func_base <- glm(outcome ~ diab + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
run_lrtest <- lrtest(glm_func, glm_func_base)
run_lrtest # same answer!



# how about 3df???
# unified model
# G is the outcome, I should simply be able to do disease*exposure

# 3df on 13:47191972
str_subset(names(figi), "6.20688121")
env<-'diab'


glm_func <- lm(X6.20688121.T.A_dose ~ outcome*diab +age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi )
glm_func_base <- lm(X6.20688121.T.A_dose ~ age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi )
run_lrtest <- lrtest(glm_func, glm_func_base)
run_lrtest # not bad. 



stargazer(run_lrtest, type = 'text', digits = 8)

out <- capture.output(stargazer(run_lrtest, 
                                title = paste0(snp_clean, " x T2D with additional adjustment covariates"), 
                                align = T, 
                                type = 'html', 
                                column.sep.width = '10pt',
                                digits = NA))

cat(paste(out, collapse = "\n"), "\n", file = paste0("~/Dropbox/test.html"), append = F)

# the results don't look as good on stargazer so... I'll just reoprt the raw thnig
# now just have to loop over the SNPs










#--------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------#

# let's compare a few SNPs
snp_list <- c("X1.71040166.G.T_dose",   "X3.185510884.A.C_dose" , "X6.20688121.T.A_dose" ,  "X7.28189411.T.C_dose",   "X8.118185025.G.A_dose" ,"X10.114754071.T.C_dose", "X10.114784926.C.T_dose", "X12.4384844.T.G_dose" ,  "X16.53811788.A.G_dose",  "X20.6442961.G.A_dose")
snp_list_clean <- gsub("_dose", "", snp_list) %>% 
            gsub("X", "", .) %>% 
            gsub("\\.", ":", .)
snp_list_clean


diab_gxescan <- readRDS("~/data/results/diab/processed/FIGI_v2.3_gxeset_diab_basic_covars_gxescan_results.rds") %>% 
  filter(SNP %in% snp_list_clean)

names(diab_gxescan)
diab_gxescan_edit <- diab_gxescan %>% 
  # dplyr::select(SNP, Subjects, Cases, chiSqG, chiSqGxE, chiSqGE, chiSq2df) %>% 
  dplyr::select(SNP, Subjects, Cases, chiSqG, chiSqGxE, chiSqGE, chiSq2df) %>% 
  mutate(chiSq3df = chiSq2df + chiSqGE,
         pval_gxescan_2df = pchisq(chiSq2df, df = 2, lower.tail = F),
         pval_gxescan_3df = pchisq(chiSq3df, df = 3, lower.tail = F))




# now calculate GLM p values, and compare (using same covariates, obviously)

env <- 'diab'
vcfid_list <- x <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", env, "_basic_covars_glm.rds"))[,'vcfid']
geno <- readRDS(paste0("/media/work/tmp/posthoc/gwis_sig_results_output_", env, ".rds")) %>% 
  dplyr::select(-contains('p0'), 
                -contains('p1'),
                -contains('p2'))
figi <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_results_exposure_main_effects.rds") %>% 
  rename_all(tolower) %>%
  rename(EUR_subset = eur_subset) %>% 
  inner_join(geno, 'vcfid') %>% 
  filter(vcfid %in% vcfid_list) %>% 
  mutate(redmeatqc2 = factor(redmeatqc2),
         vigmodlns = as.numeric(vigmodlns))


covars_1 <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))


## 2df
glm_wrapper_2df <- function(snp) {
  glm_formula <- paste0("outcome ~ ", snp, " * ", env, " + ", paste0(covars_1, collapse = " + "))
  glm_formula_base <- paste0("outcome ~ ", env, " + ", paste0(covars_1, collapse = " + "))
  
  glm_1 <- glm(glm_formula, data = figi, family = 'binomial')
  glm_2 <- glm(glm_formula_base, data = figi, family = 'binomial')
  run_lrtest <- lrtest(glm_1, glm_2)
  run_lrtest
  return(run_lrtest[, 'Pr(>Chisq)'][2])
}

out_2df <- lapply(snp_list, glm_wrapper_2df)
out_2df_final <- data.frame(SNP = snp_list_clean, pval_glm_2df = unlist(out_2df))

## 3df
glm_wrapper_3df <- function(snp) {
  lm_formula_base <- paste0(snp, " ~ outcome * ", env, " + ", paste0(covars_1, collapse = " + "))
  lm_formula <- paste0(snp, " ~ ", paste0(covars_1, collapse = " + "))
  
  lm_base <- lm(lm_formula_base, data = figi)
  lm_interaction <- lm(lm_formula, data = figi)
  run_lrtest <- lrtest(lm_base, lm_interaction)
  return(run_lrtest[, 'Pr(>Chisq)'][2])
}

out_3df <- lapply(snp_list, glm_wrapper_3df)
out_3df_final <- data.frame(SNP = snp_list_clean, pval_glm_3df = unlist(out_3df))

final_results <- inner_join(diab_gxescan_edit, out_2df_final, by = 'SNP') %>% 
  inner_join(out_3df_final, 'SNP')

saveRDS(final_results, file = "/media/work/tmp/comparison_gxescan_vs_glm_2df_3df.rds")



# get a table of chisquare p values


## 2df
glm_wrapper_2df <- function(snp) {
  glm_formula <- paste0("outcome ~ ", snp, " * ", env, " + ", paste0(covars_1, collapse = " + "))
  glm_formula_base <- paste0("outcome ~ ", env, " + ", paste0(covars_1, collapse = " + "))
  
  glm_1 <- glm(glm_formula, data = figi, family = 'binomial')
  glm_2 <- glm(glm_formula_base, data = figi, family = 'binomial')
  run_lrtest <- lrtest(glm_1, glm_2)
  run_lrtest
  return(run_lrtest[, 'Chisq'][2])
}

out_2df <- lapply(snp_list, glm_wrapper_2df)
out_2df_final <- data.frame(SNP = snp_list_clean, chiSq_glm_2df = unlist(out_2df))

## 3df
glm_wrapper_3df <- function(snp) {
  lm_formula_base <- paste0(snp, " ~ outcome * ", env, " + ", paste0(covars_1, collapse = " + "))
  lm_formula <- paste0(snp, " ~ ", paste0(covars_1, collapse = " + "))
  
  lm_base <- lm(lm_formula_base, data = figi)
  lm_interaction <- lm(lm_formula, data = figi)
  run_lrtest <- lrtest(lm_base, lm_interaction)
  return(run_lrtest[, 'Chisq'][2])
}

out_3df <- lapply(snp_list, glm_wrapper_3df)
out_3df_final <- data.frame(SNP = snp_list_clean, chiSq_glm_3df = unlist(out_3df))


final_results_chi <- inner_join(diab_gxescan_edit, out_2df_final, by = 'SNP') %>% 
  inner_join(out_3df_final, 'SNP')


saveRDS(final_results_chi, file = "/media/work/tmp/comparison_gxescan_vs_glm_2df_3df_chisq.rds")

















#-----------------------------------------------------------------------#
#just do it like they do it on the discovery channel i mean gxeswcanR
#-----------------------------------------------------------------------#


snp_list <- c("X1.71040166.G.T_dose",   "X3.185510884.A.C_dose" , "X6.20688121.T.A_dose" ,  "X7.28189411.T.C_dose",   "X8.118185025.G.A_dose" ,"X10.114754071.T.C_dose", "X10.114784926.C.T_dose", "X12.4384844.T.G_dose" ,  "X16.53811788.A.G_dose",  "X20.6442961.G.A_dose")

snp <- 'X1.71040166.G.T_dose'

glm_1 <- glm(glm_formula, data = figi, family = 'binomial')
glm_2 <- glm(glm_formula_base, data = figi, family = 'binomial')

lrtest(glm_1, glm_2)

test <- figi %>% 
  mutate(diab = as.numeric(diab))

ge_formula <- paste0(env, " ~ ", snp, " + ", paste0(covars_1, collapse = " + "))
ge_formula_base <- paste0(env, " ~ ", paste0(covars_1, collapse = " + "))
ge <- lm(ge_formula, data = test)
ge_base <- lm(ge_formula_base, data = test)

lrtest(ge, ge_base)



# ------------------------ loop it over all the stupid thing s------------------- #

## 2df
glm_wrapper_2df <- function(snp) {
  glm_formula <- paste0("outcome ~ ", snp, " * ", env, " + ", paste0(covars_1, collapse = " + "))
  glm_formula_base <- paste0("outcome ~ ", env, " + ", paste0(covars_1, collapse = " + "))
  
  glm_1 <- glm(glm_formula, data = figi, family = 'binomial')
  glm_2 <- glm(glm_formula_base, data = figi, family = 'binomial')
  run_lrtest <- lrtest(glm_1, glm_2)
  run_lrtest
  return(run_lrtest[, 'Chisq'][2])
}

out_2df <- lapply(snp_list, glm_wrapper_2df)
out_2df_final <- data.frame(SNP = snp_list_clean, chiSq_glm_2df = unlist(out_2df))

test <- figi %>% 
  mutate(diab = as.numeric(diab))


## 3df
glm_wrapper_ge <- function(snp) {
  # lm_formula_base <- paste0(snp, " ~ outcome * ", env, " + ", paste0(covars_1, collapse = " + "))
  # lm_formula <- paste0(snp, " ~ ", paste0(covars_1, collapse = " + "))
  
  ge_formula <- paste0(env, " ~ ", snp, " + ", paste0(covars_1, collapse = " + "))
  ge_formula_base <- paste0(env, " ~ ", paste0(covars_1, collapse = " + "))
  
  # lm_base <- lm(lm_formula_base, data = figi)
  # lm_interaction <- lm(lm_formula, data = figi)
  # run_lrtest <- lrtest(lm_base, lm_interaction)
  # return(run_lrtest[, 'Chisq'][2])
  
  ge_formula <- paste0(env, " ~ ", snp, " + ", paste0(covars_1, collapse = " + "))
  ge_formula_base <- paste0(env, " ~ ", paste0(covars_1, collapse = " + "))
  ge <- lm(ge_formula, data = test)
  ge_base <- lm(ge_formula_base, data = test)
  run_lrtest <- lrtest(ge, ge_base)
  return(run_lrtest[,'Chisq'][2])
}

out_ge <- lapply(snp_list, glm_wrapper_ge)
out_ge_final <- data.frame(SNP = snp_list_clean, chiSq_glm_ge= unlist(out_ge))


final_results_chi_alt <- inner_join(diab_gxescan_edit, out_2df_final, by = 'SNP') %>% 
  inner_join(out_ge_final, 'SNP') %>% 
  mutate(chiSq_glm_3df = chiSq_glm_2df + chiSq_glm_ge,
         pval_glm_alt_3df = pchisq(chiSq_glm_3df, df = 3, lower.tail = F))


saveRDS(final_results_chi_alt, file = "/media/work/tmp/comparison_gxescan_vs_glm_2df_3df_chisq_alt.rds")



