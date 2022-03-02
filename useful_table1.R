keep <- input_data$vcfid
gxe_table1 <-  readRDS(paste0("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_table1.rds"))
gxe_table1_subset <- gxe_table1 %>%
  filter(vcfid %in% keep) %>% 
  dplyr::mutate(outcome = as.factor(outcome),
                redmeatqc2 = as.factor(redmeatqc2), 
                procmeatqc2 = as.factor(procmeatqc2), 
                diab = as.factor(diab)) %>% 
  filter(sex == 0)






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
table_covariates = c("age_ref_imp", "sex", "asp_ref", "heightcm", "bmi", "hrt_ref_pm")

table1(as.formula(paste0("~ ", paste(table_covariates, collapse = "+"), "| outcome_table1")),
       data=gxe_table1_subset,
       render.continuous=my.render.cont,
       render.categorical=my.render.cat,
       render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)

