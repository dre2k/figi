
# Table for smoking subset - maybe I should include all?


# input data
ever <- readRDS(glue("/media/work/gwis_test/smk_ever/data/FIGI_v2.3_gxeset_smk_ever_basic_covars_glm.rds")) %>% 
  dplyr::pull(vcfid)

aveday <- readRDS(glue("/media/work/gwis_test/smk_aveday/data/FIGI_v2.3_gxeset_smk_aveday_basic_covars_glm.rds")) %>% 
  dplyr::pull(vcfid)

pkyr <- readRDS(glue("/media/work/gwis_test/smk_pkyr/data/FIGI_v2.3_gxeset_smk_pkyr_basic_covars_glm.rds")) %>% 
  dplyr::pull(vcfid)

xx1 <- aveday[!aveday %in% ever]

keep <- unique(c(ever, aveday, pkyr))


gxe_table1 <-  readRDS(paste0("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_table1.rds"))


gxe_table1_subset <- gxe_table1 %>%
  filter(vcfid %in% keep) %>% 
  dplyr::mutate(outcome = as.factor(outcome),
                sex = ifelse(sex == 0, "Female", "Male")) %>% 
  dplyr::mutate(across(.cols = c('fiber', 'redmeat', 'fruit', 'vegetable', "folate_tot", "calcium_tot"), ~ as.numeric(.x)))


label(gxe_table1_subset$alcoholc) = "Alcohol consumption"
label(gxe_table1_subset$age_ref_imp) = "Age (mean imputed)"
label(gxe_table1_subset$sex) = "Sex"
label(gxe_table1_subset$bmi) = "BMI (kg/m2)"
label(gxe_table1_subset$energytot_imp) = "Total energy intake (mean imputed Kcal/day)"
label(gxe_table1_subset$famhx1) = "Family history of colorectal cancer"
label(gxe_table1_subset$educ) = "Education (highest completed)"
label(gxe_table1_subset$methrswklns) = "Physical activity (MET-hr/week)"
label(gxe_table1_subset$redmeat) = "Total dietary red meat intake (servings/day)"
label(gxe_table1_subset$fruit) = "Total dietary fruit intake (servings/day)"
label(gxe_table1_subset$fiber) = "Total fiber intake (g/day)"
label(gxe_table1_subset$vegetable) = "Total dietary vegetable intake (servings/day)"
label(gxe_table1_subset$hrt_ref_pm) = "Hormone replacement therapy use"
label(gxe_table1_subset$folate_tot) = "Total folate intake (mcg/day)"
label(gxe_table1_subset$calcium_tot) = "Total calcium intake (mg/day)"





# continuous variables
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=3),
       c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

# categorical variables
my.render.cat <- function(x) {
  c("", sapply(stats.default(x),
               function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

# insert p values. hacky, requires an outcome factor with third level for p value.
# Also need to specify values by name (outcome_table1) in the 'rndr' function
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- gxe_table1_subset[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ gxe_table1_subset$outcome_table1)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(gxe_table1_subset$outcome_table1)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

# include all covariates in every descriptive statistics table
table_covariates = c("smk_ever", "age_ref_imp", "sex", "famhx1", "heightcm", "bmi", "diab", "educ", "hrt_ref_pm2", "alcoholc", "energytot_imp", "folate_tot", "calcium_tot", "redmeat", "fiber", "fruit", "vegetable")

# table_covariates = c("alcoholc", "age_ref_imp", "sex", "famhx1", "educ", "energytot_imp", "bmi", "diab",  "smk_ever", "methrswklns", "redmeat", "fruit", "vegetableqc2", "hrt_ref_pm")



table1(as.formula(paste0("~ ", paste(table_covariates, collapse = "+"), "| outcome_table1")),
       data=gxe_table1_subset,
       render.continuous=my.render.cont,
       render.categorical=my.render.cat,
       render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)



table_covariates = c("outcome", "age_ref_imp", "sex", "famhx1", "heightcm", "bmi", "diab", "educ", "hrt_ref_pm2", "alcoholc", "energytot_imp", "folate_tot", "calcium_tot", "redmeat", "fiber", "fruit", "vegetable")

table1(as.formula(paste0("~ ", paste(table_covariates, collapse = "+"), "| smk_ever_table1")),
       data=gxe_table1_subset,
       render.continuous=my.render.cont,
       render.categorical=my.render.cat,
       render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)








#-----------------------------------------------------------------------------#
# random check
#-----------------------------------------------------------------------------#


# input data
ever <- readRDS(glue("/media/work/gwis_test/smk_ever/data/FIGI_v2.3_gxeset_smk_ever_basic_covars_glm.rds")) %>% 
  dplyr::pull(vcfid)

aveday <- readRDS(glue("/media/work/gwis_test/smk_aveday/data/FIGI_v2.3_gxeset_smk_aveday_basic_covars_glm.rds")) %>% 
  dplyr::pull(vcfid)

pkyr <- readRDS(glue("/media/work/gwis_test/smk_pkyr/data/FIGI_v2.3_gxeset_smk_pkyr_basic_covars_glm.rds")) %>% 
  dplyr::pull(vcfid)

xx1 <- aveday[!aveday %in% ever]

total <- unique(c(ever, aveday, pkyr))


input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")) 

# holup - if you have aveday/pkyr you SHOULD HAVE SMK_EVER!
xx1 <- aveday[!aveday %in% ever]
wtf <- input_data %>% 
  filter(vcfid %in% xx1)
# the issue is that ATBC only had smokers, so they're included in smoking intensity analyses




x <- dplyr::filter(input_data, vcfid %in% ever) %>% 
  bind_rows(wtf)


table(x$outcome)

x_case <- dplyr::filter(x, outcome == 1)
table(x_case$smk_ever)

x_control <- dplyr::filter(x, outcome == 0)
table(x_control$smk_ever)





# input data
esubset <- readRDS(glue("/media/work/gwis_test/smk_ever/data/FIGI_v2.3_gxeset_smk_ever_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% esubset)




