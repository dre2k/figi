#=============================================================================#
# last minute stuff for manuscript
# when you get time,  you should organize code into function because you'll
# need to run this for other expsures anyway
#=============================================================================#



# nsaid vs aspirin subsets ------------------------------------------------
asp_ref <- readRDS("/media/work/gwis_test/asp_ref/data/FIGI_v2.4_gxeset_asp_ref_basic_covars_glm.rds") %>% pull(vcfid)
aspirin <- readRDS("/media/work/gwis_test/aspirin/data/FIGI_v2.4_gxeset_aspirin_basic_covars_glm.rds") %>% pull(vcfid)
length(unique(c(asp_ref, aspirin))) # all aspirin captured by asp_ref


#-----------------------------------------#
# table 1
#-----------------------------------------#
library(table1)


gxe_table1 <-  readRDS(paste0(path, "../data/FIGI_", hrc_version, "_gxeset_analysis_data_table1.rds"))

gxe_table1_subset <- gxe_table1 %>%
  filter(vcfid %in% esubset) %>% 
  dplyr::mutate(outcome = as.factor(outcome), 
                redmeatqc2 = factor(redmeatqc2, labels = c("Q1", "Q2", "Q3", "Q4")))

# appropriate labels
label(gxe_table1_subset$age_ref_imp) = "Age (mean imputed)"
label(gxe_table1_subset$sex) = "Sex"
label(gxe_table1_subset$asp_ref) = "Regular NSAIDs use"
label(gxe_table1_subset$bmi) = "BMI"
label(gxe_table1_subset$energytot_imp) = "Total energy intake (mean imputed)"
label(gxe_table1_subset$famhx) = "Family history of CRC"
label(gxe_table1_subset$educ) = "Education (highest level completed)"
label(gxe_table1_subset$smk_ever) = "Ever smoking"
label(gxe_table1_subset$diab) = "Type 2 diabetes (ever diagnosed)"
label(gxe_table1_subset$redmeatqc2) = "Total dietary red meat intake"





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
table_covariates = c("age_ref_imp", "sex", "asp_ref", "bmi", "energytot_imp", "energytot_imp", "famhx1", "educ", "smk_ever", "diab", "redmeatqc2")

table1(as.formula(paste0("~ ", paste(table_covariates, collapse = "+"), "| outcome_table1")),
       data=gxe_table1_subset,
       render.continuous=my.render.cont,
       render.categorical=my.render.cat,
       render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)



#-----------------------------------------#
# forest plot by study (not study_gxe)
#-----------------------------------------#





#-----------------------------------------#
# run analysis stratified by 
# sex + tumor site
#-----------------------------------------#

# ONLY FOR CHROMOSOME 5 finding
# get GxE model, then stratify by sex + tumor site 












