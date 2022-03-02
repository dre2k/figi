# --------------------------------------------------------------------------- #
# 5/28/2021
# additional tables for Kristina/Anna
# --------------------------------------------------------------------------- #


# ---- table1 package render functions ---- #
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
    y <- gxe_table1[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ gxe_table1$outcome_table1)$p.value
    } else {
      p <- chisq.test(table(y, droplevels(gxe_table1$outcome_table1)))$p.value
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





# ---- input data ---- #
# get IDs of individuals involved in either alcoholc_moderate OR alcoholc_heavy_vs_moderate

alcoholc_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_moderate/data/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)
alcoholc_heavy_vs_moderate_vcfid <- readRDS("/media/work/gwis_test/alcoholc_heavy_vs_moderate/data/FIGI_v2.3_gxeset_alcoholc_heavy_vs_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)

# N = 74099
keep <- unique(c(alcoholc_moderate_vcfid, alcoholc_heavy_vs_moderate_vcfid))

gxe_table1 <-  readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_table1.rds")) %>% 
  filter(vcfid %in% keep) %>% 
  mutate(outcome_table1 = fct_relevel(outcome_table1, "Cases"), 
         sex = case_when(sex == 0 ~ "Female", 
                         sex == 1 ~ "Male", 
                         TRUE ~ ""), 
         cancer_site_sum2 = factor(cancer_site_sum2), 
         alcoholc = factor(alcoholc, labels = c("Light-to-moderate drinkers (>1-28 g/d)", 
                                               "Non-drinkers (\u2264 1 g/day)", 
                                               "Heavy drinkers (>28 g/d)")), 
         methrswklns = as.numeric(methrswklns), 
         redmeatqc2 = factor(redmeatqc2, labels = c("Q1", "Q2", "Q3", "Q4")),
         fruitqc2 = factor(fruitqc2, labels = c("Q1", "Q2", "Q3", "Q4")),
         vegetableqc2 = factor(vegetableqc2, labels = c("Q1", "Q2", "Q3", "Q4")))
  

label(gxe_table1$alcoholc) = "Alcohol consumption"
label(gxe_table1$age_ref_imp) = "Age (mean imputed)"
label(gxe_table1$sex) = "Sex"
label(gxe_table1$bmi) = "BMI"
label(gxe_table1$energytot_imp) = "Total energy intake (mean imputed)"
label(gxe_table1$famhx1) = "Family history of colorectal cancer"
label(gxe_table1$educ) = "Education (highest completed)"
label(gxe_table1$smk_ever) = "Smoking"
label(gxe_table1$methrswklns) = "Physical activity (MET-hr/week)"
label(gxe_table1$redmeatqc2) = "Total dietary red meat intake"
label(gxe_table1$fruitqc2) = "Total dietary fruit intake"
label(gxe_table1$vegetableqc2) = "Total dietary vegetable intake"
label(gxe_table1$hrt_ref_pm) = "Hormone replacement therapy use"

# include all covariates in every descriptive statistics table
table_covariates = c("age_ref_imp", "sex", "asp_ref", "heightcm", "bmi", "energytot", "energytot_imp", "famhx1", "educ", "smk_ever", "hrt_ref_pm2", "diab", "calcium_totqc2", "folate_totqc2", "fiberqc2", "redmeatqc2", "procmeatqc2", "fruitqc2", "vegetableqc2", "p_diet_std")

table_covariates = c("alcoholc", "cancer_site_sum2", "age_ref_imp", "sex", "famhx1", "educ", "energytot_imp", "bmi", "diab",  "smk_ever", "methrswklns", "redmeatqc2", "fruitqc2", "vegetableqc2", "hrt_ref_pm")

uy

table1(as.formula(paste0("~ ", paste(table_covariates, collapse = "+"), "| outcome_table1")),
       data=gxe_table1,
       render.continuous=my.render.cont,
       render.categorical=my.render.cat,
       render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)






 # quick check for studies included in hRT
 hrt <- readRDS("/media/work/gwis_test/hrt_ref_pm2/data/FIGI_v2.3_gxeset_hrt_ref_pm2_basic_covars_glm.rds") %>% 
   pull(vcfid)
                
                
# input data

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% hrt)



gwas <- readRDS("/media/work/gwis_test/data/FIGI_v2.3_GWAS.rds")
table(gwas$filename)
newfoundland <- filter(gwas, filename == "newfoundland_omniquad")
                









# ----------------------------------------------- #
# calculating correlations etc for manuscript
# ----------------------------------------------- #





# requests for alcoholc x g manuscript

x <- readRDS('/media/work/gwis_test/alcoholc/output/posthoc/dosage_chr10_101351704.qs')
gwas <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gwasset_basic_covars_gxescan.rds")



gwas <- readRDS("~/data/FIGI_EpiData_rdata/FIGI_v2.3_GWAS.rds")


# GxE set
# main effect rs2300985     chr10_101476905_G_A


alcoholc_moderate <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_alcoholc_moderate_basic_covars_glm.rds") %>% 
  pull(vcfid)

tmp <- qread("/media/work/gwis_test/alcoholc/output/posthoc/dosage_chr10_101476905.qs") %>% 
  dplyr::filter(vcfid %in% alcoholc_moderate) %>% 
  inner_join(gwas, 'vcfid')

model_g <- glm(outcome ~ chr10_101476905_G_A_dose + sex + age_ref_imp + PC1 + PC2 + PC3, data = tmp, family = 'binomial')
summary(model_g)


# correlation?

# in gwas subset:
tmp <- qread("/media/work/gwis_test/alcoholc/output/posthoc/dosage_chr10_101476905.qs") %>%
  # dplyr::filter(vcfid %in% alcoholc_moderate) %>% 
  inner_join(gwas, 'vcfid')

tmp2 <- qread("/media/work/gwis_test/alcoholc/output/posthoc/dosage_chr10_101351704.qs") %>% 
  inner_join(tmp, 'vcfid')

cor(tmp2$chr10_101351704_A_G_dose, tmp2$chr10_101476905_G_A_dose)


model_g <- glm(outcome ~ chr10_101476905_G_A_dose + sex + age_ref_imp + PC1 + PC2 + PC3, data = tmp, family = 'binomial')
summary(model_g)


# correlation? 

tmp2 <- qread("/media/work/gwis_test/alcoholc/output/posthoc/dosage_chr10_101351704.qs") %>% 
  inner_join(tmp, 'vcfid')







cor(tmp2$chr10_101351704_A_G_dose, tmp2$chr10_101476905_G_A_dose)


                