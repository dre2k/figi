# Spend some time fixing your shit
# p-values are wrong for the continuous stratified odds ratio table

# so the above is fixed, now incorporate dosage based stratified table
# (no need to counts, just estimates)

rm(list = ls())

# input variables
exposure = 'folate_sup_yn'
hrc_version = 'v3.0'
path = glue("/media/work/gwis_test/{exposure}/")
covariates <- sort(c('age_ref_imp', 'sex', 'energytot_imp', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% esubset) %>% 
  mutate(fruit5qcm = as.numeric(fruit5qcm))

# inputs
data_epi = input_data
exposure
snp = '3:12041456:T:A'
snp = '6:23445253:T:A'
hrc_version
covariates
dosage = F
path = "/media/work/gwis_test/folate_sup_yn/output"
flip_allele = F




# let me  make sure this example matches RERI .. 

# input variables
exposure = 'smk_ever'
hrc_version = 'v2.3'
path = glue("/media/work/gwis_test/{exposure}/")
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))

# input data
esubset <- readRDS(glue("{path}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% pull(vcfid)
input_data <- readRDS(glue("/media/work/gwis/data/FIGI_EpiData/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid %in% esubset) %>% 
  mutate(fruit5qcm = as.numeric(fruit5qcm))

# inputs
data_epi = input_data
exposure
snp = '3:12041456:T:A'
snp = '1:151700227:G:A'
hrc_version
covariates
dosage = F
path = "/media/work/gwis_test/smk_ever/output"
flip_allele = F



# ===================================================== #
# ===================================================== #
# ===================================================== #
# stratified odds ratios ---- 
# ===================================================== #
# ===================================================== #

fit_stratified_or <- function(data_epi, exposure, snp, hrc_version, covariates, dosage = F, path, flip_allele = F) {
  
  mod = glue_collapse(covariates, sep = "+")
  wdir = glue("{path}/posthoc")
  
  # SNP information
  snp_info <- unlist(strsplit(snp, split = ":"))
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  snpfix <- snpname_clean(snp)
  snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))
  
  
  
  
  # read dosage file
  data_dose <- qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')
  
  if (flip_allele == T) {
    # flip dosages
    # data[, snp] <- abs(data[, paste0(snp)] - 2)
    data[, snpfix] <- abs(data[, snpfix] - 2)
    # flip genotype probabilities
    pp <- data[,paste0(snpfix_short, "_p2")]
    data[,paste0(snpfix_short, "_p2")] <- data[, paste0(snpfix_short, "_p0")]
    data[,paste0(snpfix_short, "_p0")] <- pp
    # assign ref/alt allele
    ref_allele <- snp_info[4]
    alt_allele <- snp_info[3]
  } else {
    ref_allele <- snp_info[3]
    alt_allele <- snp_info[4]
  }
  
  
  
  # create data subset
  tmp1 = data[, c('outcome', exposure, covariates)]
  tmp2 = data[, grepl(paste(paste0(snpfix_short, c("_dose", "_p0", "_p1", "_p2")), collapse = "|"), names(data))]
  names(tmp2) <- c("dosage", "p0", "p1", "p2")
  ds_tmp = cbind(tmp1, tmp2) %>%
    na.omit(.[, c(exposure, 'outcome', covariates, 'p1', 'p2', 'dosage')]) # from data.table. remove rows with NA in these columns)
  
  
  
  
  #---- Calculate counts for each cell ----------------
  Ncaco = data.frame(snps=snpfix_short,matrix(NA,ncol=6*2))
  colnames(Ncaco) = c('snps',paste0('Co',1:2),paste0('Ca',1:2),
                      paste0('Co1',1:2),paste0('Ca1',1:2),
                      paste0('Co2',1:2),paste0('Ca2',1:2))
  rownames(Ncaco) = snpfix_short
  
  Ncaco[snpfix_short,c(paste0('Co',1:2),paste0('Ca',1:2))] = t(table(ds_tmp[,c('outcome',exposure)]))
  Ncaco[snpfix_short,c(paste0('Co1',1:2),paste0('Ca1',1:2))] = c(t(tapply(ds_tmp$p1, ds_tmp[,c('outcome',exposure)], sum, na.rm=T)))
  Ncaco[snpfix_short,c(paste0('Co2',1:2),paste0('Ca2',1:2))] = c(t(tapply(ds_tmp$p2, ds_tmp[,c('outcome',exposure)], sum, na.rm=T)))
  
  
  
  
  #== calculate counts for G=0 and put counts into format ca/co
  for(i in 1:2){
    eval(parse(text=paste0('Ncaco$caco0',i,"=paste0(round(Ncaco$Ca",i,'-Ncaco$Ca1',i,'-Ncaco$Ca2',i,
                           ",1),'/',round(Ncaco$Co",i,'-Ncaco$Co1',i,'-Ncaco$Co2',i,',1))')))
    eval(parse(text=paste0('Ncaco$caco1',i,'=paste0(round(Ncaco$Ca1',i,",1),'/',round(Ncaco$Co1",i,",1))")))
    eval(parse(text=paste0('Ncaco$caco2',i,'=paste0(round(Ncaco$Ca2',i,",1),'/',round(Ncaco$Co2",i,",1))")))
  }
  
  
  
  
  
  #-- Fit unrestricted model --------
  # ds_tmp[,exposure] = factor(ds_tmp[,exposure])
  tmp = as.vector(Model.all.new(ds_tmp, mod, exposure))
  tmp = as.vector(Model.all.new.dosage(ds_tmp, mod, exposure))
  
  res.pool.un = res.pool.g = res.pool.e = NULL
  res.pool.un = rbind(res.pool.un,data.frame(snpfix_short,exposure,t(tmp$GE)))
  res.pool.e = rbind(res.pool.e,data.frame(snpfix_short,exposure,t(tmp$E)))
  res.pool.g = rbind(res.pool.g,data.frame(snpfix_short,exposure,t(tmp$G)))

  
  ## organize results ######
  ## need to automate elvl eventually.. 
  if(is.factor(ds_tmp[,exposure])) {
    elvl <- c(0:(nlevels(ds_tmp[, exposure])-1))
  } else {
    elvl <- c(0,1)
  }
  
  # always set to genotype (have different function for per allelic dosage.. )
  glvl <- c(0:2)
  glvl <- c(0:1)
  
  ##== GE stratified results ########
  colnames(res.pool.un) = c('snp','env',
                            paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
                            'beta0p1','se0p1',
                            paste0('beta',elvl[-1],'p1'),paste0('se',elvl[-1],'p1'))
  res.pool.un = format_res(res.pool.un)
  
  ##== stratified by G results #######
  colnames(res.pool.g) = c('snp','env',paste0(rep(c('beta0','se0','beta1','se1'),each=1),
                                              rep(elvl[-1],4)))
  res.pool.g = format_res(res.pool.g)
  
  ##== stratified by E results #######
  colnames(res.pool.e) = c('snp','env',paste0('beta1',elvl),paste0('se1',elvl))
  res.pool.e = format_res(res.pool.e)
  
  
  ##== Put into table
  OR.tab = ORtab(snpfix_short,
                 elvl=elvl,
                 glvl=glvl,
                 res=res.pool.un)
  pval.tab = ptab(snpfix_short,elvl=elvl,glvl=glvl,res=res.pool.un)
  
  ORg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('OR0',elvl[-1])]),
                                         as.character(res.pool.g[,paste0('OR1',elvl[-1])])))), ncol=2)
  pg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('Ppval0',elvl[-1])]),
                                        as.character(res.pool.g[,paste0('Ppval1',elvl[-1])])))), ncol=2)
  colnames(ORg.tab) = colnames(pg.tab) = c(paste0('p',glvl))
  
  ORe.tab = matrix(as.character(unlist(c(res.pool.e[,paste0('OR1',elvl)]))), ncol=2)
  pe.tab  = matrix(as.character(unlist(c(res.pool.e[,paste0('Ppval1',elvl)]))),ncol=2)
  
  
  
  #=== write into table
  if(is.factor(ds_tmp[, exposure])) {
    est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
                      ORg.tab[1,],pg.tab[1,]),
                rbind(ORe.tab[1,1],pe.tab[1,1],ORe.tab[1,2],pe.tab[1,2],rep(NA,1),rep(NA,1)),
                N0 = c(Ncaco[snpfix_short,'caco01'],NA,Ncaco[snpfix_short,'caco02'],rep(NA,3)),
                N1 = c(Ncaco[snpfix_short,'caco11'],NA,Ncaco[snpfix_short,'caco12'],rep(NA,3)),
                N2 = c(Ncaco[snpfix_short,'caco21'],NA,Ncaco[snpfix_short,'caco22'],rep(NA,3)))
    
    colnames(est) <- c(paste0(ref_allele, ref_allele),
                       glue("{ref_allele}/{alt_allele}, {alt_allele}/{alt_allele}"),
                       glue("{ref_allele}/{alt_allele}, {alt_allele}/{alt_allele}"),
                       "N0", "N1", "N2")
    exposure_level <- levels(ds_tmp[,exposure])
    rownames(est) <- c(paste0(exposure,"=",exposure_level[1]),
                       "p0",
                       paste0(exposure,"=",exposure_level[2]),
                       "p1",
                       exposure,
                       "p")
    
  } else {
    est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
                      ORg.tab[1,],pg.tab[1,]),
                rbind(ORe.tab[1,],pe.tab[1,],ORe.tab[2,],pe.tab[2,],rep(NA,2),rep(NA,2)))
    colnames(est) <- c(paste0(ref_allele, ref_allele),
                       paste0(ref_allele, alt_allele),
                       paste0(alt_allele, alt_allele),
                       paste0(ref_allele, alt_allele),
                       paste0(alt_allele, alt_allele))
    exposure_level <- levels(ds_tmp[,exposure])
    rownames(est) <- c(paste0(exposure,"=Ref"),
                       "p0",
                       paste0(exposure,"=UnitInc"),
                       "p1",
                       exposure,
                       "p")
  }
  
  saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  
  
  # # slightly modify output in case anyone wants per allele effects rather than p1/p2
  # if(dosage == T) {
  #   
  #   tmp = as.vector(Model.all.new.dosage(ds_tmp,mod,exposure))
  #   res.pool.un = res.pool.g = res.pool.e = NULL
  #   res.pool.un = rbind(res.pool.un,data.frame(snpfix_short,exposure,t(tmp$GE)))
  #   
  #   elvl=c(0:1) ; glvl = c(0,1)
  #   tmp2 = NULL
  #   tmp2 = rbind(tmp2,data.frame(snpfix_short,exposure,t(tmp$GE)))
  #   
  #   colnames(tmp2) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
  #                      'beta0p1','beta0p2','se0p1','se0p2')
  #   tmp2 = format_res(tmp2)
  #   
  #   res.pool.e = rbind(res.pool.e,data.frame(snpfix_short,exposure,t(tmp$E)))
  #   colnames(res.pool.e) = c('snp',
  #                            'env',
  #                            paste0('beta1',elvl),paste0('se1',elvl))
  #   tmp2 <- format_res(res.pool.e)
  #   
  #   
  #   ORe.tab = matrix(as.character(unlist(c(tmp2[,paste0('OR1',elvl)]))),ncol=2)
  #   pe.tab  = matrix(as.character(unlist(c(tmp2[,paste0('Ppval1',elvl)]))),ncol=2,byrow=T)
  #   est2 <- rbind(ORe.tab[1,1],pe.tab[1,1],ORe.tab[1,2],pe.tab[1,2], NA, NA)
  #   
  #   final_out <- est %>%
  #     dplyr::select(-c(4,5)) %>%
  #     add_column(est2, .after = 3)
  #   
  #   colnames(final_out) <- c(paste0(ref_allele, ref_allele),
  #                            paste0(ref_allele, alt_allele),
  #                            paste0(alt_allele, alt_allele),
  #                            paste0(alt_allele, " Allelic Dosage"),
  #                            "N0", "N1", "N2")
  #   # this is only for FACTOR variables (might have to modify when you run Q4 variables)
  #   exposure_level <- levels(data[,exposure])
  #   rownames(final_out) <- c(paste0(exposure,"=",exposure_level[1]),
  #                            "p0",
  #                            paste0(exposure,"=",exposure_level[2]),
  #                            "p1",
  #                            exposure,
  #                            "p")
  #   saveRDS(final_out, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}_dosage.rds"))
  #   return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}_dosage.rds"))
  # } else {
  #   colnames(est) <- c(paste0(ref_allele, ref_allele),
  #                      paste0(ref_allele, alt_allele),
  #                      paste0(alt_allele, alt_allele),
  #                      paste0(ref_allele, alt_allele),
  #                      paste0(alt_allele, alt_allele),
  #                      "N0", "N1", "N2")
  #   exposure_level <- levels(data[,exposure])
  #   rownames(est) <- c(paste0(exposure,"=",exposure_level[1]),
  #                      "p0",
  #                      paste0(exposure,"=",exposure_level[2]),
  #                      "p1",
  #                      exposure,
  #                      "p")
  #   
  #   saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  #   return(glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
  #   # saveRDS(est, file = glue("{wdir}/stratified_oddsratio_{exposure}_{hrc_version}_{snpfix_short}.rds"))
  # }
  
}






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
model <- glm(glue("outcome ~ {exposure}*{snp_new} + {glue_collapse(covariates, sep = '+')}"), 
             family = binomial(link = "logit"), data = data)
coef_names <- names(coef(model))
coef_exposure <- grep(exposure, coef_names)[1]
coef_snp <- grep(snp_new, coef_names)[1]
coef_interaction <- grep(exposure, coef_names)[2]
reri_est = epi.interaction(model = model, coef = c(coef_exposure, 
                                                   coef_snp, coef_interaction), param = "product", type = "RERI", 
                           conf.level = 0.95)
coef_keep <- coef_names[c(coef_exposure, coef_snp, coef_interaction)]
cov.mat <- vcov(model)
V2 = cov.mat[coef_keep, coef_keep]
reri_se = deltamethod(~exp(x1 + x2 + x3) - exp(x1) - exp(x2) + 
                        1, mean = c(coef(model)[coef_exposure], coef(model)[coef_snp], 
                                    coef(model)[coef_interaction]), cov = V2)
reri_pval = format.pval(2 * pnorm(-abs(reri_est[1, 1]/reri_se)), 
                        digits = 4)
value = interactionR(model, exposure_names = c(exposure, 
                                               snp_new), ci.type = "delta", ci.level = 0.95, em = F, 
                     recode = F)
out <- interactionR_table2(value, pvalue = reri_pval)
saveRDS(out, file = glue("{wdir}/reri_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
}


# ========================== #
# interaction plot ----
# ========================== #

iplot_wrapper <- function(data_epi, exposure, hrc_version, snp, covariates, path, flip_allele = F) {
  
  wdir = glue("{path}/posthoc")
  
  # SNP - need to change the SNP name to match the file (because of binarydosage output)
  tmp <- unlist(strsplit(snp, ":"))
  chr <- tmp[1]
  bp <- as.numeric(tmp[2])
  ref <- tmp[3]
  alt <- tmp[4]
  
  # convert 12:50610976:C:T to chr12_50610976_C_T_dose
  snpname_clean <- function(x) {
    tmp <- gsub("\\:", "\\_", x)
    # tmp <- gsub("X", "chr", tmp)
    tmp <- glue("chr{tmp}_dose")
    return(tmp)
  }
  
  snpfix <- snpname_clean(snp)
  
  #  need to merge with dosage information
  data_dose <- qread(glue("{wdir}/dosage_chr{chr}_{bp}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')
  
  if (flip_allele == T) {
    snp_old <- snpfix
    snp_tmp <- unlist(strsplit(snpfix, split = "_"))
    chr <- snp_tmp[1]
    bp <- snp_tmp[2]
    a1 <- snp_tmp[3]
    a2 <- snp_tmp[4]
    snp_new <- glue("{chr}_{bp}_{a2}_{a1}_dose_flipped")
    data[[snp_new]] <- abs(2-data[, snp_old])
  } else {
    snp_new <- snpfix
  }
  
  model <- glm(glue("outcome ~ {exposure}*{snp_new} + {glue_collapse(covariates, sep = '+')}"), family = binomial(link = "logit"), data = data)
  
  png(glue("{wdir}/interaction_plot_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}.png"), height = 720, width = 1280)
  if (is.factor(data[,exposure])) {
    print(interact_plot(model, modx = !! exposure , pred = !! snp_new, plot.points = F, interval = T, outcome.scale = 'link', y.label = 'predicted log odds') + theme(text = element_text(size = 26)))
    # johnson_neyman(model, pred = folate_totqc2, modx = chr2_55255870_C_T, alpha = 0.05)
  } else {
    print(interact_plot(model, modx = !! exposure , pred = !! snp_new, plot.points = F, interval = T, modx.values = c(0,1), outcome.scale = 'link', y.label = 'predicted log odds') + theme(text = element_text(size = 26)))
  }
  dev.off()
  
  return(glue("{wdir}/interaction_plot_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}.png"))
  # saveRDS(out, file = glue("{output_dir}reri_{exposure}_{snp}_{glue_collapse(sort(covariates), sep = '_')}.rds"))
}
















# ======================================= #
# Gxe stratified by alcoholc ----
# ======================================= #

data_epi <- input_data
exposure 
snp = '3:12041456:T:A'
snp = '6:23445253:T:A'
covariates
strata = 'alcoholc'
method = 'chiSqGxE'
flip_allele = F
path = '/media/work/gwis_test/folate_sup_yn/output/'


function (data_epi, exposure, snp, covariates, strata = c("sex", 
                                                          "study_design", "cancer_site_sum2"), method = c("chiSqGxE", 
                                                                                                          "two-step", "chiSqCase", "chiSq2df", "chiSq3df"), flip_allele = F, 
          path) 
{
  
  
  
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
  
  if (flip_allele == T) {
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
  
  
  strata <- match.arg(strata)
  method <- match.arg(method)
  
  covariates_nostrata <- paste0(covariates[!covariates %in% strata], collapse = " + ")
  
  data[, "strata_num"] <- as.numeric(factor(data[, strata])) - 1
  data[, "exposure_num"] = as.numeric(data[, exposure])
  
  out <- list()
  out_all <- fit_gxe(data, exposure, snp_new, covariates)
  out[["all"]] <- out_all
  number_of_levels <- nlevels(factor(data[, strata]))
  
  
  for (level in seq(number_of_levels) - 1) {
    if (strata == "cancer_site_sum2") {
      index_vector <- which(data[, "strata_num"] == level | 
                              data[, "outcome"] == 0)
    }
    else {
      index_vector <- which(data[, "strata_num"] == level)
    }
    out_level <- fit_gxe(data[index_vector, ], exposure, 
                         snp_new, covariates_nostrata)
    out[[paste0(strata, "_", as.character(level))]] <- out_level
  }
  
  
  
  list_of_glms <- lapply(out, function(x) x[[1]])
  list_of_samplesizes <- lapply(list_of_glms, function(x) paste0(c("Ca=", 
                                                                   "Co="), rev(as.character(table(x$model$outcome))), collapse = ","))
  coefs <- lapply(list_of_glms, function(x) (exp(coef(x))))
  if (strata == "sex") {
    col_label = paste0(c("All", "Female", "Male"), " (", 
                       list_of_samplesizes, ")")
  }
  else if (strata == "study_design") {
    col_label = paste0(c("All", "Cohort", "Case-Control"), 
                       " (", list_of_samplesizes, ")")
  }
  else if (strata == "cancer_site_sum2") {
    col_label = paste0(c("All", "Proximal", "Distal", "Rectal"), 
                       " (", list_of_samplesizes, ")")
  } 
  else if (strata == "alcoholc") {
    col_label = paste0(c("All", "1-28g/day", "nondrinker", ">28g/day"), 
                       " (", list_of_samplesizes, ")")
  }
  
  
  
  list_of_chisq <- lapply(out, function(x) x[[2]])
  if (method %in% c("chiSqGxE", "two-step", "chiSqCase")) {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[1]], 
                                                                               df = 1, lower.tail = F), format = "e", digits = 4)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("GxE term LRtest p = ", gxe_pvalues))
  }
  else if (method == "chiSq2df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[2]], 
                                                                               df = 2, lower.tail = F), format = "e", digits = 4)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("2DF LRtest p = ", gxe_pvalues))
  }
  else if (method == "chiSq3df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq((x[[2]] + 
                                                                                  x[[3]]), df = 3, lower.tail = F), format = "e", digits = 4)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("3DF LRtest p = ", gxe_pvalues))
  }
  out_html <- stargazer_helper(list_of_glms, title = paste0(gsub("\\_", "\\\\_", strata), " stratified ", gsub("\\_", "\\\\_", snp_new), " x ", gsub("\\_", "\\\\_", exposure)), column.labels = col_label, coef = coefs, notes = notes, single.row = T)
  
  cat(paste(out_html, collapse = "\n"), "\n", file = glue("~/Dropbox/gxe_models_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}_stratified_by_{strata}.html"), append = F)
  return(glue("~/Dropbox/gxe_models_{exposure}_{hrc_version}_{snpfix}_{glue_collapse(sort(covariates), sep = '_')}_stratified_by_{strata}.html"))
}




