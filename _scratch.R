
# trying to further adjusted 2-step plot ----------------------------------
# e.g. multiple methods = splitting overall alpha accordingly

# (use diab as example)

data = input_data
exposure
path = path
hrc_version = 'v2.3'
simplem_step1_statistic = 'chiSqG'
covars = covariates
binsToPlot = 7
stats_step1 = 'chiSqG'
sizeBin0 = 5
alpha = 0.05
m = 1000000
meff_method = 'gao'
gao_pca_cutoff = 0.995
number_of_snps






simplem_wrap <- function(data, exposure, hrc_version, covariates, simplem_step1_statistic, path, meff_method, gwas_snps = NULL, gao_pca_cutoff = 0.995) {
  
  # create list of bin SNP dosage data.frames
  # FYI - 'glue' doesn't seem regex friendly
  files_input <- mixedsort(list.files(glue("{path}/data/expectation_bin_dosages/"), pattern = paste0(paste0("expectation_based_snplist_", exposure, "_", hrc_version, "_", simplem_step1_statistic, "_bin"), "(?:.+)", "output.rds"), full.names = T)) 
  files_list <- map(files_input, ~ readRDS(.x)) 
  
  # ------------------------------------------- #
  # function to filter GWAS SNPs if desired
  remove_snps <- function(zz, gwas_snps_vector) {
    zznames <- substr(names(zz), 1, nchar(names(zz)) - 4)
    zz_index <- !zznames %in% gwas_snps_vector
    zz_out <- zz[, zz_index]
    return(zz_out)
  }
  
  # filter out GWAS SNPs if provided
  # be careful - make sure gwas_snps is provided as chr:bp - so i can change it to Xchr.bp
  if (!is.null(gwas_snps)) {
    # GWAS SNPs to remove (vector of chr:bp)
    #gwas_snps <- readRDS("../data/conditioning_snps_v20200930_filter_GWAS_SNPS_1000k.rds")
    gwas_snps_v2 <- paste0("X", gsub(":", ".", gwas_snps))
    
    # remove from bins
    files_list <- map(files_list, ~ remove_snps(.x, gwas_snps_v2))
    
    # remove from data
    data2 <- data %>%
      dplyr::filter(!SNP2 %in% gwas_snps)
  } else {
    data2 <- data
  }    
  
  # output correlation matrix for each bin 
  # overall, no longer doing it by chromosome
  snps_cor <- map(files_list, ~ cor(.x[,-1])) # '-1' to remove vcfid column
  
  # ------------------------------------------ #
  number_of_snps <- map_int(files_list, ~ ncol(.x[,-1]))
  
  # calculate effective number of tests, output results in vector
  if(meff_method == "gao") {
    number_of_tests <- map_dbl(snps_cor, ~ poolr::meff(R = .x, method = meff_method, C = gao_pca_cutoff))
  } else if(meff_method == "lea") {
    number_of_tests <- map_dbl(snps_cor, ~ meff_lea(s = .x))
  } else {
    number_of_tests <- map_dbl(snps_cor, ~ poolr::meff(R = .x, method = meff_method))
  }
  
  # helper for writing filename
  if(meff_method == "gao") {
    meff_method_out = glue("gao_{gao_pca_cutoff}")
  } else {
    meff_method_out = meff_method
  }
  
  # helper for writing 'gwas' if you filtered gwas 
  if(!is.null(gwas_snps)) {
    meff_method_out = glue(meff_method_out, "_excludeGWAS")
  }
  






# folate total included ---------------------------------------------------

x1 <- readRDS("/media/work/gwis_test/folate_totqc2/data/FIGI_v3.0_gxeset_folate_totqc2_basic_covars_glm.rds") %>% pull(vcfid)
x2 <- readRDS("/media/work/gwis_test/folate_dietqc2/data/FIGI_v3.0_gxeset_folate_dietqc2_basic_covars_glm.rds") %>% pull(vcfid)
x3 <- readRDS("/media/work/gwis_test/folate_sup_yn/data/FIGI_v3.0_gxeset_folate_sup_yn_basic_covars_glm.rds") %>% pull(vcfid)


x <- unique(c(x1, x2, x3))

input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_v3.0_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% x) 

table(input_data$outcome)


# interaction plot (for diab paper i believe) -----------------------------


snp <- c("8:118185025:G:A")
dos <- qread("/media/work/gwis_test/diab/output/posthoc/dosage_chr8_118185025.qs")
out <- inner_join(input_data, dos, 'vcfid') %>% 
  mutate(diab_num = as.numeric(diab)) %>% 
  filter(!study_gxe %in% c("PPS3", "PPS4"))


data_epi = out
function (data_epi, exposure, snp, covariates, path) 
{
  wdir <- glue("{path}/posthoc")
  snp_info <- unlist(strsplit(snp, split = ":"))
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  data_dose <- qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, "vcfid")
  model_check <- glm(glue("outcome ~ {exposure}*{snpfix} + {glue_collapse(covariates, sep = '+')}"), 
                     family = "binomial", data = data)
  
  data <- out
  if (model_check[[1]][snpfix] < 0) {
    snp_old <- snpfix
    snp_tmp <- unlist(strsplit(snpfix, split = "_"))
    chr <- snp_tmp[1]
    bp <- snp_tmp[2]
    a1 <- snp_tmp[3]
    a2 <- snp_tmp[4]
    snp_new <- glue("{chr}_{bp}_{a2}_{a1}_dose_flipped")
    data[[snp_new]] <- abs(2 - data[, snp_old])
  }
  else {
    snp_new <- snpfix
  }
  
  # data = out
  
  model <- glm(glue("outcome ~ {exposure}*{snp_new} + {glue_collapse(covariates, sep = '+')}"), 
               family = binomial(link = "logit"), data = data)
  coef_names <- names(coef(model))
  coef_exposure <- grep(exposure, coef_names)[1]
  coef_snp <- grep(snp_new, coef_names)[1]
  coef_interaction <- grep(exposure, coef_names)[2]
  reri_est = epi.interaction(model = model, coef = c(coef_exposure, 
                                                     coef_snp, coef_interaction), param = "product",
                             conf.level = 0.95)
  coef_keep <- coef_names[c(coef_exposure, coef_snp, coef_interaction)]
  cov.mat <- vcov(model)
  V2 = cov.mat[coef_keep, coef_keep]
  reri_se = deltamethod(~exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 
                          1, mean = c(coef(model)[coef_exposure], coef(model)[coef_snp], 
                                      coef(model)[coef_interaction]), cov = V2)
  reri_pval = format.pval(2 * pnorm(-abs(reri_est[1, 1]/reri_se)), 
                          digits = 4)
  reri_pval = 0
  value = interactionR(model, exposure_names = c(exposure, 
                                                 snp_new), ci.type = "delta", ci.level = 0.95, em = F, 
                       recode = F)
  out2 <- interactionR_table2(value, pvalue = reri_pval)
  saveRDS(out, file = glue("{wdir}/reri_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
}




# stratified odds ratio tables for specific use case ----------------------
# I used this to create stratified odds ratio for diabetes paper (it was a mess, but got it done)
# had to do it this way because in the results I included per allele effect, and to get pvalue i had to do it manually (using SE to calculate wald test)

dos <- qread("/media/work/gwis_test/diab/output/posthoc/dosage_chr8_118185025.qs")
out <- inner_join(input_data, dos, 'vcfid') %>% 
  mutate(diab_num = as.numeric(diab)) %>% 
  filter(!study_gxe %in% c("PPS3", "PPS4"))


data_epi <- out
exposure
# snp = "13:47191972:G:A"
covariates



function (data_epi, exposure, snp, hrc_version, covariates, dosage = F, 
          path, flip_allele = F) 
{
  
  mod = glue_collapse(covariates, sep = "+")
  wdir = glue("{path}/posthoc")
  snp_info <- unlist(strsplit(snp, split = ":"))
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  data_dose <- qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, "vcfid")
  data <- out
  if (flip_allele == T) {
    data[, snpfix] <- abs(data[, snpfix] - 2)
    pp <- data[, paste0(snpfix_short, "_p2")]
    data[, paste0(snpfix_short, "_p2")] <- data[, paste0(snpfix_short, 
                                                         "_p0")]
    data[, paste0(snpfix_short, "_p0")] <- pp
    ref_allele <- snp_info[4]
    alt_allele <- snp_info[3]
  }
  else {
    ref_allele <- snp_info[3]
    alt_allele <- snp_info[4]
  }
  tmp1 = data[, c("outcome", exposure, covariates)]
  tmp2 = data[, grepl(paste(paste0(snpfix_short, c("_dose", 
                                                   "_p0", "_p1", "_p2")), collapse = "|"), names(data))]
  names(tmp2) <- c("dosage", "p0", "p1", "p2")
  ds_tmp = cbind(tmp1, tmp2) %>% na.omit(.[, c(exposure, "outcome", 
                                               covariates, "p1", "p2", "dosage")])
  Ncaco = data.frame(snps = snpfix_short, matrix(NA, ncol = 6 * 
                                                   2))
  colnames(Ncaco) = c("snps", paste0("Co", 1:2), paste0("Ca", 
                                                        1:2), paste0("Co1", 1:2), paste0("Ca1", 1:2), paste0("Co2", 
                                                                                                             1:2), paste0("Ca2", 1:2))
  rownames(Ncaco) = snpfix_short
  Ncaco[snpfix_short, c(paste0("Co", 1:2), paste0("Ca", 1:2))] = t(table(ds_tmp[, 
                                                                                c("outcome", exposure)]))
  Ncaco[snpfix_short, c(paste0("Co1", 1:2), paste0("Ca1", 1:2))] = c(t(tapply(ds_tmp$p1, 
                                                                              ds_tmp[, c("outcome", exposure)], sum, na.rm = T)))
  Ncaco[snpfix_short, c(paste0("Co2", 1:2), paste0("Ca2", 1:2))] = c(t(tapply(ds_tmp$p2, 
                                                                              ds_tmp[, c("outcome", exposure)], sum, na.rm = T)))
  for (i in 1:2) {
    eval(parse(text = paste0("Ncaco$caco0", i, "=paste0(round(Ncaco$Ca", 
                             i, "-Ncaco$Ca1", i, "-Ncaco$Ca2", i, ",1),'/',round(Ncaco$Co", 
                             i, "-Ncaco$Co1", i, "-Ncaco$Co2", i, ",1))")))
    eval(parse(text = paste0("Ncaco$caco1", i, "=paste0(round(Ncaco$Ca1", 
                             i, ",1),'/',round(Ncaco$Co1", i, ",1))")))
    eval(parse(text = paste0("Ncaco$caco2", i, "=paste0(round(Ncaco$Ca2", 
                             i, ",1),'/',round(Ncaco$Co2", i, ",1))")))
  }
  tmp = as.vector(Model.all.new.dosage(ds_tmp, mod, exposure))
  
  
  
  
  
  res.pool.un = res.pool.g = res.pool.e = NULL
  res.pool.un = rbind(res.pool.un, data.frame(snpfix_short, 
                                              exposure, t(tmp$GE)))
  res.pool.e = rbind(res.pool.e, data.frame(snpfix_short, exposure, 
                                            t(tmp$E)))
  res.pool.g = rbind(res.pool.g, data.frame(snpfix_short, exposure, 
                                            t(tmp$G)))
  if (is.factor(ds_tmp[, exposure])) {
    elvl <- c(0:(nlevels(ds_tmp[, exposure]) - 1))
  }
  else {
    elvl <- c(0, 1)
  }
  glvl <- c(0:2)
  colnames(res.pool.un) = c("snp", "env", paste0("beta", elvl[-1], 
                                                 "p0"), paste0("se", elvl[-1], "p0"), "beta0p1", "beta0p2", 
                            "se0p1", "se0p2", paste0("beta", elvl[-1], "p1"), paste0("se", 
                                                                                     elvl[-1], "p1"), paste0("beta", elvl[-1], "p2"), 
                            paste0("se", elvl[-1], "p2"))
  res.pool.un = format_res(res.pool.un)
  colnames(res.pool.g) = c("snp", "env", paste0(rep(c("beta0", 
                                                      "se0", "beta1", "se1", "beta2", "se2"), each = 1), rep(elvl[-1], 
                                                                                                             6)))
  res.pool.g = format_res(res.pool.g)
  colnames(res.pool.e) = c("snp", "env", paste0("beta1", elvl), 
                           paste0("se1", elvl), paste0("beta2", elvl), paste0("se2", 
                                                                              elvl))
  res.pool.e = format_res(res.pool.e)
  OR.tab = ORtab(snpfix_short, elvl = elvl, glvl = glvl, res = res.pool.un)
  pval.tab = ptab(snpfix_short, elvl = elvl, glvl = glvl, res = res.pool.un)
  ORg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[, 
                                                                 paste0("OR0", elvl[-1])]), as.character(res.pool.g[, 
                                                                                                                    paste0("OR1", elvl[-1])]), as.character(res.pool.g[, 
                                                                                                                                                                       paste0("OR2", elvl[-1])])))), ncol = 3)
  pg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[, 
                                                                paste0("Ppval0", elvl[-1])]), as.character(res.pool.g[, 
                                                                                                                      paste0("Ppval1", elvl[-1])]), as.character(res.pool.g[, 
                                                                                                                                                                            paste0("Ppval2", elvl[-1])])))), ncol = 3)
  colnames(ORg.tab) = colnames(pg.tab) = c(paste0("p", glvl))
  ORe.tab = matrix(as.character(unlist(c(res.pool.e[, paste0("OR1", 
                                                             elvl)], res.pool.e[, paste0("OR2", elvl)]))), ncol = 2)
  pe.tab = matrix(as.character(unlist(c(res.pool.e[, paste0("Ppval1", 
                                                            elvl)], res.pool.e[, paste0("Ppval2", elvl)]))), ncol = 2)
  if (is.factor(ds_tmp[, exposure])) {
    est = cbind(rbind(OR.tab[1, ], pval.tab[1, ], OR.tab[2, 
    ], pval.tab[2, ], ORg.tab[1, ], pg.tab[1, ]), rbind(ORe.tab[1, 
    ], pe.tab[1, ], ORe.tab[2, ], pe.tab[2, ], rep(NA, 
                                                   2), rep(NA, 2)), N0 = c(Ncaco[snpfix_short, "caco01"], 
                                                                           NA, Ncaco[snpfix_short, "caco02"], rep(NA, 3)), N1 = c(Ncaco[snpfix_short, 
                                                                                                                                        "caco11"], NA, Ncaco[snpfix_short, "caco12"], rep(NA, 
                                                                                                                                                                                          3)), N2 = c(Ncaco[snpfix_short, "caco21"], NA, Ncaco[snpfix_short, 
                                                                                                                                                                                                                                               "caco22"], rep(NA, 3)))
    colnames(est) <- c(paste0(ref_allele, ref_allele), paste0(ref_allele, 
                                                              alt_allele), paste0(alt_allele, alt_allele), paste0(ref_allele, 
                                                                                                                  alt_allele), paste0(alt_allele, alt_allele), "N0", 
                       "N1", "N2")
    exposure_level <- levels(ds_tmp[, exposure])
    rownames(est) <- c(paste0(exposure, "=", exposure_level[1]), 
                       "p0", paste0(exposure, "=", exposure_level[2]), "p1", 
                       exposure, "p")
  }
  else {
    est = cbind(rbind(OR.tab[1, ], pval.tab[1, ], OR.tab[2, 
    ], pval.tab[2, ], ORg.tab[1, ], pg.tab[1, ]), rbind(ORe.tab[1, 
    ], pe.tab[1, ], ORe.tab[2, ], pe.tab[2, ], rep(NA, 
                                                   2), rep(NA, 2)))
    colnames(est) <- c(paste0(ref_allele, ref_allele), paste0(ref_allele, 
                                                              alt_allele), paste0(alt_allele, alt_allele), paste0(ref_allele, 
                                                                                                                  alt_allele), paste0(alt_allele, alt_allele))
    exposure_level <- levels(ds_tmp[, exposure])
    rownames(est) <- c(paste0(exposure, "=Ref"), "p0", paste0(exposure, 
                                                              "=UnitInc"), "p1", exposure, "p")
  }
  saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
}





# removing PPS from T2D paper ---------------------------------------------
snps <- c("13:47191972:G:A", "8:118185025:G:A")

## 13:47191972:G:A (see table 1 in Niki's manuscript)
## original p-values match, go ahead and filter PPS3/4
dos <- qread("/media/work/gwis_test/diab/output/posthoc/dosage_chr13_47191972.qs")
out <- inner_join(input_data, dos, 'vcfid') %>% 
  mutate(diab_num = as.numeric(diab)) 

model_base <- glm(outcome ~ diab + chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_gxe <- glm(outcome ~ diab * chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_gwas <- glm(outcome ~ diab + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out, family = 'binomial')
model_base_eg <- lm(diab_num ~ age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out)
model_eg <- lm(diab_num ~ chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out)

lrtest(model_base, model_gwas)
lrtest(model_base_eg, model_eg)
lrtest(model_base, model_gxe)
lrtest(model_gwas, model_gxe)

model_gxe <- glm(outcome ~ diab * chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi, data = out, family = 'binomial')
model_gwas <- glm(outcome ~ diab + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe + bmi, data = out, family = 'binomial')
lrtest(model_gwas, model_gxe)




model_base <- glm(outcome ~ diab*chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out , family = 'binomial')
summary(model_base)

exp(7.568e-02)
exp(7.568e-02 -2.374e-01)




model_base <- glm(outcome ~ diab*chr13_47191972_G_A_dose + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = out , family = 'binomial', 
                  contrasts = list(diab = 1))

print(high_fat, X = TRUE)







# gt package explore (nice) -----------------------------------------------


input_data %>%
  # filter(study %in% studies) %>%
  dplyr::select(study_gxe, outcome) %>% 
  group_by(outcome, study_gxe) %>% 
  dplyr::summarise(n = n()) %>% 
  pivot_wider(names_from = outcome , values_from = n) %>% 
  gt()


input_data %>%
  filter(grepl("CCFR", study_gxe)) %>%
  dplyr::select(study_gxe, study_site, outcome) %>% 
  group_by(outcome, study_gxe, study_site) %>% 
  dplyr::summarise(n = n())

table(input_data$study_gxe)






# forest plot -------------------------------------------------------------

data_epi = input_data
exposure
covariates = c("age_ref_imp", "sex")
hrc_version
path
forest_height = 17
forest_width = 8.5
strata = 'all'
categorical = F





tmp_fp <- function (data_epi, exposure, covariates, hrc_version, path, 
          forest_height = 17, forest_width = 8.5, funnel_height = 8, 
          funnel_width = 8.5, strata = "all", categorical = T) 
{
  covariates <- covariates[which(covariates != "study_gxe")]
  if (strata == "female") {
    data_epi <- dplyr::filter(data_epi, sex == 0)
    covariates <- covariates[which(!covariates %in% c("sex"))]
  }
  else if (strata == "male") {
    data_epi <- dplyr::filter(data_epi, sex == 1)
    covariates <- covariates[which(!covariates %in% c("sex"))]
  }
  else if (strata == "proximal") {
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 
                                "proximal" | data_epi$outcome == 0)
  }
  else if (strata == "distal") {
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 
                                "distal" | data_epi$outcome == 0)
  }
  else if (strata == "rectal") {
    data_epi <- dplyr::filter(data_epi, cancer_site_sum2 == 
                                "rectal" | data_epi$outcome == 0)
  }
  
  
  if (categorical == T) {
    drops <- data.frame(table(data_epi$study_gxe, data_epi$outcome, 
                              data_epi[, exposure])) %>% dplyr::filter(Freq == 
                                                                         0) %>% dplyr::pull(Var1) %>% unique(.)
    data_epi <- data_epi %>% dplyr::filter(!study_gxe %in% 
                                             drops)
  }
  else if (categorical == F) {
    drops <- data.frame(table(data_epi$study_gxe, data_epi$outcome)) %>% 
      dplyr::filter(Freq <= 5) %>% dplyr::pull(Var1) %>% 
      unique(.)
    data_epi <- data_epi %>% dplyr::filter(!study_gxe %in% 
                                             drops)
  }
  
  
  
  
  model_formula <- glue("outcome ~ {exposure} + {glue_collapse(sort(covariates), sep = '+')}")
  
  study_design <- data_epi %>% dplyr::select(study_gxe, study_design) %>% 
    filter(!duplicated(.)) %>% arrange(study_gxe)
  
  glm_out <- data_epi %>% 
    tidyr::nest(data = -study_gxe) %>% 
    dplyr::mutate(fit = purrr::map(data, ~glm(model_formula, data = .x, family = "binomial")), 
                  tidied = purrr::map(fit, ~tidy(.x)), 
                  quality = purrr::map(fit, ~glance(.x))) %>% 
    dplyr::select(-data, -fit) %>% tidyr::unnest(tidied) %>% 
    tidyr::unnest(quality) %>% dplyr::filter(grepl(exposure, 
                                                   term)) %>% dplyr::arrange(study_gxe)
  meta_input <- get_counts_outcome_by_group(data_epi, "outcome", 
                                            "study_gxe") %>% mutate(study_gxe = as.character(study_gxe)) %>% 
    inner_join(study_design, "study_gxe") %>% full_join(glm_out, 
                                                        "study_gxe")
  results_meta <- meta::metagen(estimate, std.error, data = meta_input, 
                                studlab = paste(study_gxe), fixed = TRUE, random = TRUE, 
                                method.tau = "SJ", hakn = TRUE, prediction = TRUE, sm = "OR", 
                                subgroup = study_design)
  plot_title = glue(strata, " (N=", sum(glm_out$nobs), ")", 
                    "\n{model_formula}")
  png(glue("{path}/forest_plot_{exposure}_{hrc_version}_", 
           glue_collapse(covariates, "_"), "_{strata}.png"), height = forest_height, 
      width = forest_width, units = "in", res = 150)
  meta::forest(results_meta, layout = "JAMA", leftcols = c("studlab", 
                                                           "Control", "Case", "N", "effect", "ci", "w.random"), 
               digits.addcols = 0, study.results = T, prediction = F, 
               col.random = "red", xlim = c(0.25, 4))
  grid.text(plot_title, 0.5, 0.98, gp = gpar(cex = 1))
  dev.off()
  png(glue("{path}/funnel_plot_{exposure}_{hrc_version}_", 
           glue_collapse(covariates, "_"), "_{strata}.png"), height = funnel_height, 
      width = funnel_width, units = "in", res = 150)
  meta::funnel(results_meta, sm = "OR", studlab = T, pos = 4, 
               col.random = "red")
  dev.off()
}








# forest plot ver 2 -------------------------------------------------------------
# some Q4 exposures issue where some studies have 2 exposure levels instead of 4
# don't filter empty cells that's all 

# have to run this bit by bit, not working currently
data_epi = input_data
exposure
covariates = c("age_ref_imp", "sex", "energytot_imp")
hrc_version
path = glue("/media/work/gwis_test/{exposure}/output/posthoc/")
forest_height = 17
forest_width = 9
strata = 'all'
categorical = F

funnel_height = 8
funnel_width = 8.5


tmp_fp2 <- function (data_epi,
                    exposure,
                    covariates,
                    hrc_version,
                    path,
                    forest_height = 17,
                    forest_width = 8.5,
                    funnel_height = 8,
                    funnel_width = 8.5,
                    strata = "all",
                    categorical = T)
{
  
  
  covariates <- covariates[which(covariates != "study_gxe")]
  if (strata == "female") {
    data_epi <- dplyr::filter(data_epi, sex == 0)
    covariates <- covariates[which(!covariates %in% c("sex"))]
  }
  else if (strata == "male") {
    data_epi <- dplyr::filter(data_epi, sex == 1)
    covariates <- covariates[which(!covariates %in% c("sex"))]
  }
  else if (strata == "proximal") {
    data_epi <- dplyr::filter(data_epi,
                              cancer_site_sum2 ==
                                "proximal" | data_epi$outcome == 0)
  }
  else if (strata == "distal") {
    data_epi <- dplyr::filter(data_epi,
                              cancer_site_sum2 ==
                                "distal" | data_epi$outcome == 0)
  }
  else if (strata == "rectal") {
    data_epi <- dplyr::filter(data_epi,
                              cancer_site_sum2 ==
                                "rectal" | data_epi$outcome == 0)
  }
  # if (categorical == T) {
  #   drops <- data.frame(table(data_epi$study_gxe, data_epi$outcome,
  #                             data_epi[, exposure])) %>% dplyr::filter(Freq ==0) %>% 
  #     dplyr::pull(Var1) %>% unique(.)
  #   data_epi <- data_epi %>% dplyr::filter(!study_gxe %in%
  #                                            drops)
  # }
  # else if (categorical == F) {
  #   drops <-
  #     data.frame(table(data_epi$study_gxe, data_epi$outcome)) %>%
  #     dplyr::filter(Freq <= 5) %>% dplyr::pull(Var1) %>%
  #     unique(.)
  #   data_epi <- data_epi %>% dplyr::filter(!study_gxe %in%
  #                                            drops)
  # }
  model_formula <-
    glue("outcome ~ {exposure} + {glue_collapse(sort(covariates), sep = '+')}")
  
  study_design <-
    data_epi %>% dplyr::select(study_gxe, study_design) %>%
    filter(!duplicated(.)) %>% arrange(study_gxe)
  
  glm_out <- data_epi %>% 
    tidyr::nest(data = -study_gxe) %>%
    dplyr::mutate(fit = purrr::map(data, ~ glm(model_formula, data = .x, family = "binomial")),
                  tidied = purrr::map(fit, ~ tidy(.x)),
                  quality = purrr::map(fit, ~ glance(.x))) %>%
    dplyr::select(-data,-fit) %>% 
    tidyr::unnest(tidied) %>%
    tidyr::unnest(quality) %>% 
    dplyr::filter(grepl(exposure, term)) %>% 
    dplyr::arrange(study_gxe)
  
  meta_input <- get_counts_outcome_by_group(data_epi, "outcome", "study_gxe") %>% 
    mutate(study_gxe = as.character(study_gxe)) %>%
    inner_join(study_design, "study_gxe") %>% 
    full_join(glm_out,"study_gxe")
  
  results_meta <-
    meta::metagen(
      estimate,
      std.error,
      data = meta_input,
      studlab = paste(study_gxe),
      fixed = TRUE,
      random = TRUE,
      method.tau = "SJ",
      hakn = TRUE,
      prediction = TRUE,
      sm = "OR",
      subgroup = study_design
    )
  
  
  plot_title = glue(strata, " (N=", sum(glm_out$nobs), ")",
                    "\n{model_formula}")
  png(
    glue(
      "{path}/forest_plot_{exposure}_{hrc_version}_",
      glue_collapse(covariates, "_"),
      "_{strata}.png"
    ),
    height = forest_height,
    width = forest_width,
    units = "in",
    res = 150
  )
  meta::forest(
    results_meta,
    layout = "JAMA",
    leftcols = c("studlab",
                 "Control", "Case", "N", "effect", "ci", "w.random"),
    digits.addcols = 0,
    study.results = T,
    prediction = F,
    col.random = "red",
    xlim = c(0.5, 2)
  )
  grid.text(plot_title, 0.5, 0.98, gp = gpar(cex = 1))
  dev.off()
  png(
    glue(
      "{path}/funnel_plot_{exposure}_{hrc_version}_",
      glue_collapse(covariates, "_"),
      "_{strata}.png"
    ),
    height = funnel_height,
    width = funnel_width,
    units = "in",
    res = 150
  )
  meta::funnel(
    results_meta,
    sm = "OR",
    studlab = T,
    pos = 4,
    col.random = "red"
  )
  dev.off()
}








# table1 code to create one off tables ------------------------------------


# use subset of indiviudals included in either OR redmeat/procmeat (i think if redemat, also in procmeat)


# input variables
# exposure = 'calcium_totqc2'
exposure = 'calcium_dietqc2'
hrc_version = 'v3.0'
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")

# good overlap, not perfect
# rm_subset <- readRDS(glue("/media/work/gwis_test/redmeatqcm_v2/data/FIGI_{hrc_version}_gxeset_redmeatqcm_v2_basic_covars_glm.rds")) %>% pull(vcfid)
# pm_subset <- readRDS(glue("/media/work/gwis_test/procmeatqcm_v2/data/FIGI_{hrc_version}_gxeset_procmeatqcm_v2_basic_covars_glm.rds")) %>% pull(vcfid)
# subset <- unique(c(rm_subset, pm_subset))

subset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)

gxe_table1_subset <-  readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_table1.rds")) %>% 
  mutate(sex = factor(sex, labels = c("Female", "Male")), 
         bmic3 = factor(bmic3, labels = c("Normal", "Obese", "Overweight"))) %>% 
  filter(vcfid%in% subset) %>%
  # mutate(redmeatqcm_v2 = as.numeric(redmeatqcm), 
  #        procmeatqcm_v2 = as.numeric(procmeatqcm))
  mutate(across(.cols = c(calcium_tot, folate_tot, fiber, fruit, vegetable, redmeat, procmeat),  ~ as.numeric(.x)),
         across(.cols = c(calcium_totqc2, folate_totqc2, fiberqc2, fruitqc2, vegetableqc2, redmeatqc2, procmeatqc2),  ~ as.factor(.x)))








# ------------- # 
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=3),
       c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x),
               function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}


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


table_covariates = c("redmeatqcm_v2", "procmeatqcm_v2", "age_ref_imp", "sex", "heightcm", "bmi", "bmic3", "energytot", "energytot_imp", "famhx1", "educ", "asp_ref", "smk_ever", "hrt_ref_pm2", "diab", "calcium_totqc2", "folate_totqc2", "fiberqc2", "redmeatqc2", "procmeatqc2", "fruitqc2", "vegetableqc2", "alcoholc", "methrswklns")

# table_covariates = c("redmeatqcm_v2", "procmeatqcm_v2", "age_ref_imp", "sex", "heightcm", "bmi", "bmic3", "energytot", "energytot_imp", "famhx1", "educ", "asp_ref", "smk_ever", "hrt_ref_pm2", "diab", "calcium_tot", "folate_tot", "fiber", "redmeat", "procmeat", "fruit", "vegetable")

table1(as.formula(paste0("~ ", paste(table_covariates, collapse = "+"), "| outcome_table1")),
       data=gxe_table1_subset,
       render.continuous=my.render.cont,
       render.categorical=my.render.cat,
       render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)







# table1 code to create one off tables - BY EXPOSURE ------------------------------------


# use subset of indiviudals included in either OR redmeat/procmeat (i think if redemat, also in procmeat)


# input variables
# exposure = 'calcium_totqc2'
exposure = 'calcium_dietqc2'
hrc_version = 'v3.0'
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")

# good overlap, not perfect
# rm_subset <- readRDS(glue("/media/work/gwis_test/redmeatqcm_v2/data/FIGI_{hrc_version}_gxeset_redmeatqcm_v2_basic_covars_glm.rds")) %>% pull(vcfid)
# pm_subset <- readRDS(glue("/media/work/gwis_test/procmeatqcm_v2/data/FIGI_{hrc_version}_gxeset_procmeatqcm_v2_basic_covars_glm.rds")) %>% pull(vcfid)
# subset <- unique(c(rm_subset, pm_subset))

subset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)

gxe_table1_subset <-  readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_table1.rds")) %>% 
  mutate(sex = factor(sex, labels = c("Female", "Male")), 
         bmic3 = factor(bmic3, labels = c("Normal", "Obese", "Overweight"))) %>% 
  filter(vcfid%in% subset) %>%
  # mutate(redmeatqcm_v2 = as.numeric(redmeatqcm), 
  #        procmeatqcm_v2 = as.numeric(procmeatqcm))
  mutate(across(.cols = c(calcium_tot, folate_tot, fiber, fruit, vegetable, redmeat, procmeat),  ~ as.numeric(.x)),
         across(.cols = c(calcium_totqc2, folate_totqc2, fiberqc2, fruitqc2, vegetableqc2, redmeatqc2, procmeatqc2),  ~ as.factor(.x)))








# ------------- # 
my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=3),
       c("", "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}

my.render.cat <- function(x) {
  c("", sapply(stats.default(x),
               function(y) with(y, sprintf("%d (%0.0f %%)", FREQ, PCT))))
}

rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- gxe_table1_subset[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y) & length(table(gxe_table1_subset[,paste0(exposure, '_table1')])) <= 2) {
      p <- t.test(y ~ gxe_table1_subset[, paste0(exposure, '_table1')])$p.value
    } else if (is.numeric(y) & length(table(gxe_table1_subset[,paste0(exposure, '_table1')])) > 2) {
      p <- anova(lm(y ~ gxe_table1_subset[, paste0(exposure, '_table1')]))$`Pr(>F)`
    } else {
      p <- chisq.test(table(y, droplevels(gxe_table1_subset[, paste0(exposure, '_table1')])))$p.value
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


table_covariates = c("outc", "redmeatqcm_v2", "procmeatqcm_v2", "age_ref_imp", "sex", "heightcm", "bmi", "bmic3", "energytot", "energytot_imp", "famhx1", "educ", "asp_ref", "smk_ever", "hrt_ref_pm2", "diab", "calcium_totqc2", "folate_totqc2", "fiberqc2", "redmeatqc2", "procmeatqc2", "fruitqc2", "vegetableqc2", "alcoholc", "methrswklns")

# table_covariates = c("redmeatqcm_v2", "procmeatqcm_v2", "age_ref_imp", "sex", "heightcm", "bmi", "bmic3", "energytot", "energytot_imp", "famhx1", "educ", "asp_ref", "smk_ever", "hrt_ref_pm2", "diab", "calcium_tot", "folate_tot", "fiber", "redmeat", "procmeat", "fruit", "vegetable")

table1(as.formula(paste0("~ ", paste(table_covariates, collapse = "+"), "| ", paste0(exposure, "_table1"))),
       data=gxe_table1_subset,
       render.continuous=my.render.cont,
       render.categorical=my.render.cat,
       render=rndr, render.strat=rndr.strat, overall = F, droplevels = F)





# investigating sex-study specific quartiles ------------------------------


table(input_data$study_gxe, input_data$outcome)

# investigate red meat ... 

mec1 <- filter(input_data, study_gxe == "MEC_1")

mec1 %>% 
  group_by(outcome) %>% 
  summarise(mm = mean(as.numeric(calcium_tot))) %>% 
  as.data.frame(.)

out <- mec1 %>% 
  group_by(sex) %>% 
  summarise(quantile = scales::percent(c(0.25, 0.5, 0.75)),
            calcium_tot = quantile(as.numeric(calcium_tot), c(0.25, 0.5, 0.75)))


cc <- filter(mec1, sex == 0 & outcome == 0) %>% pull(calcium_tot) %>% as.numeric()
ccm <- filter(mec1, sex == 1 & outcome == 0) %>% pull(calcium_tot) %>% as.numeric()

cc <- filter(mec1, sex == 0) %>% pull(calcium_tot) %>% as.numeric()
ccm <- filter(mec1, sex == 1) %>% pull(calcium_tot) %>% as.numeric()



ccq = quantile(cc, c(0.25, 0.5, 0.75), na.rm = T)
ccmq = quantile(ccm, c(0.25, 0.5, 0.75), na.rm = T)



mec1f <- filter(mec1, sex == 0) %>% 
  mutate(q4 = cut(as.numeric(calcium_tot), c(-Inf, ccq, Inf), label = c("0", "1", "2", "3")))
mec1m <- filter(mec1, sex == 1) %>% 
  mutate(q4 = cut(as.numeric(calcium_tot), c(-Inf, ccmq, Inf), label = c("0", "1", "2", "3")))

out2 <- bind_rows(mec1f, mec1m) %>% 
  filter(redmeat != )
table(out2$q4, out2$outcome)

# ------------ red meat ----------- #

# real answer

d1 <- filter(mec1, sex == 0)
table(d1$redmeatqcm)
d2 <- filter(mec1, sex == 1)
table(d2$redmeatqcm)



cc <- filter(mec1, sex == 0 & outcome == 0) %>% pull(redmeat) %>% as.numeric()
ccm <- filter(mec1, sex == 1 & outcome == 0) %>% pull(redmeat) %>% as.numeric()

cc <- filter(mec1, sex == 0) %>% pull(redmeat) %>% as.numeric()
ccm <- filter(mec1, sex == 1) %>% pull(redmeat) %>% as.numeric()

ccq = quantile(cc, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)
ccmq = quantile(ccm, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)



mec1f <- filter(mec1, sex == 0) %>% 
  mutate(q4 = cut(as.numeric(redmeat), c(-Inf, ccq, Inf), label = c("0", "1", "2", "3")))
mec1m <- filter(mec1, sex == 1) %>% 
  mutate(q4 = cut(as.numeric(redmeat), c(-Inf, ccmq, Inf), label = c("0", "1", "2", "3")))

mec1f <- filter(mec1, sex == 0) %>% 
  mutate(q4 = cut(as.numeric(redmeat), ccq, label = c("0", "1", "2", "3")))
mec1m <- filter(mec1, sex == 1) %>% 
  mutate(q4 = cut(as.numeric(redmeat), ccmq, label = c("0", "1", "2", "3")))



out2 <- bind_rows(mec1f, mec1m)
table(out2$q4, out2$outcome)


# -------------- #
# i'm beginning to think this was done by STUDY (not study/platform)

# for fucks sake, just check males ok 

mec <- filter(input_data, study == "MECC")

d2 <- filter(mec, sex == 1) %>% 
  filter(study_gxe == "MECC_1")
table(d2$redmeatqcm)



ccm <- filter(mec, sex == 1 & outcome == 0) %>% pull(redmeat) %>% as.numeric()
ccm <- pull(mec, redmeat) %>% as.numeric()
ccm <- filter(mec, sex == 1) %>% pull(redmeat) %>% as.numeric()

ccmq = quantile(ccm, c(0, 0.25, 0.5, 0.75, 1), na.rm = T)

mecm <- filter(mec, sex == 1) %>% 
  mutate(q4 = cut(as.numeric(redmeat), ccmq, label = c("0", "1", "2", "3")))

table(mecm %>% filter(study_gxe == "MECC_1") %>% pull(q4))



