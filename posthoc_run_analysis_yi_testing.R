#=============================================================================#
# FIGI GWIS Post hoc analyses
# 03/16/2020
#
# Create the following:
# - summary of significant and suggestive findings, by analysis method
# - stratified plots of predicted logits by strata (E and G)
# - stratified table of odds ratios. Binary vars only?
# - locuszoom plots of loci
# - functional annotation plots
#=============================================================================#
library(tidyverse)
library(data.table)
library(ggplot2)
library(qqman)
library(table1)
library(meta)
library(rlang)
library(broom)
library(effects)
library(figifs)
library(DT)
library(grid)
library(sjPlot)
library(stargazer)
library(forcats)
library(kableExtra)
rm(list = ls())

env = 'hrt_ref_pm2'
covs = c('age_ref_imp', 'study_gxe', 'pc1', 'pc2', 'pc3')

args <- commandArgs(trailingOnly=T)
env <- args[1] # ex: asp_ref
mod <- args[2]
covs <- c(args[3:length(args)]) # space delim list of covars, turns into vector

dir.create(paste0("/media/work/tmp/posthoc/", env), recursive = T, showWarnings = F)


#-----------------------------------------------------------------------------#
# input data ----
#-----------------------------------------------------------------------------#

# snps (dose + probability)
geno <- readRDS(paste0("/media/work/tmp/posthoc/gwis_sig_results_output_", env, ".rds"))
dat <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", env, "_basic_covars_glm.rds")) %>% 
  inner_join(geno, 'vcfid') %>% 
  rename_at(vars(c("PC1", "PC2", "PC3")), tolower) # for neatness

# flip E coding to risk increase
if(env %in% c("hrt_ref_pm2", "ep_ref_pm_gxe", "eo_ref_pm_gxe")) {
  # dat <- mutate(dat, !!env := abs(!!(sym(env)) - 1))
  dat[, env] <- abs(dat[, env] - 1)
}





#-----------------------------------------------------------------------------#
# stratified predicted logit plots ----
#-----------------------------------------------------------------------------#

# GLM wrap function. Outcome is ALWAYS ALWAYS ALWAYS 'outcome'
plot_wrapper <- function(dat, exposure, snp, covariates = c()) {
  
  # check if exposure variable is binary, fit models accordingly (categorical vs numeric)
  factor_or_numeric <- function(exposure) {
    if(length(table(exposure)) <= 4) {
      return(as.factor)
    } else {
      return(as.numeric)
    }}
  
  # change SNP name to cleaner format
  snp_rename <- (gsub('\\.', '_', snp) %>%  gsub('X', 'chr', .) %>% gsub("_dose", "", .))
  
  tmp <- dat %>% 
    mutate(!!snp_rename := as.numeric(!!sym(snp))) %>% 
    mutate_at(vars(contains(exposure)), factor_or_numeric(exposure))
  
  # fit model, output plots and odds ratios data.frame
  mod_formula <- as.formula(paste0('outcome ~ ', exposure , "*", snp_rename, "+", paste0(covariates, collapse = "+")))
  mod <- glm(mod_formula, data = tmp, family = 'binomial')
  
  mod_list <- list(x=c(0,1,2))
  names(mod_list) <- snp_rename
  model_eff <- effect(paste0(exposure,":", snp_rename),
                      mod,
                      se = TRUE,
                      confidence.level = 0.95,
                      typical = mean,
                      family = 'binomial',
                      xlevels = mod_list)
  
  # output plots as png files
  oldw <- getOption("warn")
  options(warn = -1)
  
  png(paste0("/media/work/tmp/posthoc/", exposure, "/effects_plot_", exposure, "_", snp_rename, "_", paste0(sort(covariates), collapse = "_"), "_e.png"), height = 4, width = 6, units = 'in', res = 150)
  print({
    plot(model_eff, snp_rename, multiline = T, confint = list(style = 'bands'), type = 'link', ylab = 'predicted outcome log-odds')
  })
  dev.off()
  
  png(paste0("/media/work/tmp/posthoc/", exposure, "/effects_plot_", exposure, "_", snp_rename, "_", paste0(sort(covariates), collapse = "_"), "_g.png"), height = 4, width = 6, units = 'in', res = 150)
  print({
    plot(model_eff, exposure, multiline = T, confint = list(style = 'auto'), type = 'link', ylab = 'predicted outcome log-odds')
  })
  dev.off()
  
  options(warn = oldw)
}

snp_list <- names(dat)[grepl("^X\\d{1,2}.*dose$", names(dat))]
snp_list

# run plot wrapper
sapply(snp_list, function(x) plot_wrapper(dat, exposure = exposure, snp = x, covariates = covariates))


#-----------------------------------------------------------------------------#
# Yi's script to generate stratified odds ratios
#-----------------------------------------------------------------------------#

# suppose you generated reports, and extracted dosages/genotype probabilities for significant and suggestive SNPs
# now you're interested in creating stratified odds ratios tables to better visualize the interaction
# you can use genotype probabilities to get genotype specific ORs for E and vice versa. Let's look at Yi's code

# first, organize dataset of exposure and genetic probabilities for a single SNP as an example
# let's use as an example X6.117823508.T.C



#---------------------------#
# BUNCH OF TESTING STUFF
# # input data
# geno <- readRDS(paste0("/media/work/tmp/posthoc/gwis_sig_results_output_", exposure, ".rds"))
# dat <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds")) %>% 
#   inner_join(geno, 'vcfid') %>% 
#   rename_at(vars(c("PC1", "PC2", "PC3")), tolower) # for neatness
# 
# # create vector of snp variables to apply wrapper over
# snp_list <- names(dat)[grepl("X\\d{1,2}", names(dat))] # anything starting with a X# or X##
# snp_list <- snp_list[grepl("?dose", snp_list)] # get dose only (shortcut because regex no good)
# snp_list <- gsub("_dose", "", snp_list) # take out _dose
# snp_list
# 
# sapply(snp_list, function(x) glm_wrapper(dat, exposure = exposure, snp = x, covariates = covariates))
# 
# 
# 
# 
# # using X6.117823508 as an example
# hrt <- readRDS("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_hrt_ref_pm2_basic_covars_glm.rds")
# snp_values <- readRDS("/media/work/tmp/posthoc/gwis_sig_results_output_hrt_ref_pm2.rds") 
# # dplyr::select(-X6.117823508.T.C_dose)
# # names(snp_values) <- c("vcfid", "dosage", "p0", "p1", "p2")
# 
# dat <- inner_join(hrt, snp_values, 'vcfid') %>% 
#   rename_at(vars(c("PC1", "PC2", "PC3")), tolower) 
# 
# 
# 
# 
# 
covs = c('age_ref_imp','study_gxe','pc1','pc2','pc3')
mod = 'age_ref_imp+study_gxe+pc1+pc2+pc3+study_gxe:pc1+study_gxe:pc2+study_gxe:pc3+study_gxe:age_ref_imp+study_gxe:p1 + study_gxe:p2+ study_gxe:hrt_ref_pm2'
mod = 'hrt_ref_pm2 + p1 + p2 + hrt_ref_pm2:p1 + hrt_ref_pm2:p2 + age_ref_imp + study_gxe + pc1 + pc2 + pc3'
env ='hrt_ref_pm2' ## E
snps='X6.148292041.T.C' ##  SNP
s='X6.148292041.T.C' ##  SNP


res.pool.un = res.pool.g = res.pool.e = NULL
Ncaco = data.frame(snps=snps,matrix(NA,ncol=6*2))
colnames(Ncaco) = c('snps',paste0('Co',1:2),paste0('Ca',1:2),
                    paste0('Co1',1:2),paste0('Ca1',1:2),
                    paste0('Co2',1:2),paste0('Ca2',1:2))
rownames(Ncaco) = snps
e=env


#== stratified analysis for gxe ==================
# Yi functions
Model.all.new <- function(data,mod,env)
{
  if(length(table(data$p1))<=1 | length(table(data$p2))<=1 | 
     var(2*data$p2+data$p1,na.rm=T)==0 | 
     sum((data$p2+data$p1)>1.1,na.rm=T)>0) {
    res = NA
  } else {
    data[,env]=factor(data[,env])
    mModel <- paste0('outcome~',env,'*p1+',env,'*p2+',mod)
    tmp <- summary(glm(mModel, family='binomial', data=data))
    COV <- tmp$cov.unscaled
    env.idx = env
    if(is.factor(data[,env])) env.idx = paste0(env,levels(data[,env])[-1])
    
    # stratified by GE  
    eg0 <- tmp$coef[env.idx,c(1,2)]
    beta.eg1 <- tmp$coef[env.idx,1]+tmp$coef['p1',1]+tmp$coef[paste0(env.idx,':p1'),1]
    beta.eg2 <- tmp$coef[env.idx,1]+tmp$coef['p2',1]+tmp$coef[paste0(env.idx,':p2'),1]
    e0g <- tmp$coef[c('p1','p2'),c(1,2)]
    I3 =  rep(1,3)
    se.eg1 = se.eg2 = NULL
    for(i in 1:(nlevels(data[,env])-1))
    {
      covs = COV[c(env.idx[i],'p1',paste0(env.idx[i],':p1')),c(env.idx[i],'p1',paste0(env.idx[i],':p1'))]		
      se.eg1 = c(se.eg1,sqrt(t(I3)%*%covs%*%I3))
      covs = COV[c(env.idx[i],'p2',paste0(env.idx[i],':p2')),c(env.idx[i],'p2',paste0(env.idx[i],':p2'))]		
      se.eg2 = c(se.eg2,sqrt(t(I3)%*%covs%*%I3))
    }
    GE<- c(eg0,e0g,beta.eg1,se.eg1,beta.eg2,se.eg2)
    
    # stratified by G 
    g0 <-tmp$coef[env.idx,1:2]
    if(length(env.idx)==1){
      idx.g1 = c(env.idx,paste0(env.idx,':p1'))
      idx.g2 = c(env.idx,paste0(env.idx,':p2'))
      g1 = sum(tmp$coef[idx.g1,1],na.rm=T)
      g2 = sum(tmp$coef[idx.g2,1],na.rm=T)
      Is1 =  rep(1,length(idx.g1))
      Is2 =  rep(1,length(idx.g2))
    }else{
      est.g1 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p1'),1])
      est.g2 <- data.frame(tmp$coef[env.idx,1],tmp$coef[paste0(env.idx,':p2'),1])
      g1 = rowSums(est.g1,na.rm=T)
      g2 = rowSums(est.g2,na.rm=T)
      Is1 =  rep(1,ncol(est.g1))
      Is2 =  rep(1,ncol(est.g2))
    }
    
    se.g1 = se.g2 = NULL
    for(i in 1:(nlevels(data[,env])-1))
    {
      covs1 = COV[c(env.idx[i],paste0(env.idx[i],':p1')),c(env.idx[i],paste0(env.idx[i],':p1'))]  
      covs2 = COV[c(env.idx[i],paste0(env.idx[i],':p2')),c(env.idx[i],paste0(env.idx[i],':p2'))]  
      se.g1 <- c(se.g1,sqrt(t(Is1)%*%covs1%*%Is1))
      se.g2 <- c(se.g2,sqrt(t(Is2)%*%covs2%*%Is2))
    }
    
    G<- c(g0,g1,se.g1,g2,se.g2)
    
    # stratified by E
    
    e0 <-tmp$coef[c('p1','p2'),1:2]
    if(length(env.idx)==1){
      idx.e11 = c('p1',paste0(env.idx,':p1'))
      idx.e12 = c('p2',paste0(env.idx,':p2'))
      e11 = sum(tmp$coef[idx.e11,1],na.rm=T)
      e12 = sum(tmp$coef[idx.e12,1],na.rm=T)
      Is1 =  rep(1,length(idx.e11))
      Is2 =  rep(1,length(idx.e12))
    }else{
      est.e11 <- data.frame(tmp$coef['p1',1],tmp$coef[paste0(env.idx,':p1'),1])
      est.e12 <- data.frame(tmp$coef['p2',1],tmp$coef[paste0(env.idx,':p2'),1])
      e11 = rowSums(est.e11,na.rm=T)
      e12 = rowSums(est.e12,na.rm=T)
      Is1 =  rep(1,ncol(est.e11))
      Is2 =  rep(1,ncol(est.e12))
    }
    
    se.e11 = se.e12 = NULL
    for(i in 1:(nlevels(data[,env])-1))
    {
      covs1 = COV[c('p1',paste0(env.idx[i],':p1')),c('p1',paste0(env.idx[i],':p1'))]  
      covs2 = COV[c('p2',paste0(env.idx[i],':p2')),c('p2',paste0(env.idx[i],':p2'))]  
      se.e11 <- c(se.e11,sqrt(t(Is1)%*%covs1%*%Is1))
      se.e12 <- c(se.e12,sqrt(t(Is2)%*%covs2%*%Is2))
    }
    E <- c(e0[1,1],e11,e0[1,2],se.e11,e0[2,1],e12,e0[2,2],se.e12)
  }
  res<-list(GE=GE,G=G,E=E)
}






ORtab = function(x,elvl,glvl,res) {
  ORs = data.frame(matrix(NA,nrow=length(elvl),ncol=length(glvl)))
  colnames(ORs) = paste0('p',glvl)
  rownames(ORs) = paste0('OR',elvl)
  ORs[1,1] = 1
  for(c in colnames(ORs))
    for(r in rownames(ORs))
      if(!(c=='p0' & r==paste0('OR',elvl[1]))) 
        ORs[r,c] =  as.character(res[res$snp %in% x,paste0(r,c)])
  ORs      
} 
ptab = function(x,elvl,glvl,res) {
  pval = data.frame(matrix(NA,nrow=length(elvl),ncol=length(glvl)))
  colnames(pval) = paste0('p',glvl)
  rownames(pval) = elvl
  pval[1,1] = 1
  for(c in colnames(pval))
    for(r in rownames(pval))
      if(!(c=='p0' & r==elvl[1])) 
        pval[r,c] =  as.numeric(as.vector(res[res$snp %in% x,paste0('pval',r,c)]))
  pval = apply(pval,2,function(y) paste0('P= ',formatC(y,format='g',digit=2)))    
}

format.res =function(res){
  betas = colnames(res)[grep('beta',colnames(res),fixed=T)]
  ses = colnames(res)[grep('se',colnames(res),fixed=T)]
  ORs = sapply(seq(length(betas)),function(x,betas,ses,res) 
  {
    or = paste0(rnd2(exp(res[,betas[x]])),' (',
                rnd2(exp(res[,betas[x]]-qnorm(0.975)*res[,ses[x]])),'-',
                rnd2(exp(res[,betas[x]]+qnorm(0.975)*res[,ses[x]])),')')
  },betas=betas,ses=ses,res=res)
  if(nrow(res)==1)  ORs = data.frame(t(ORs)) 
  colnames(ORs) = sub('beta','OR',betas)
  pvals = sapply(seq(length(betas)),function(x,betas,ses,res) 
  {pval = 2*pnorm(-abs(res[,betas[x]]/res[,ses[x]]))},betas=betas,ses=ses,res=res)
  if(nrow(res)==1) pvals = data.frame(t(pvals)) 
  colnames(pvals) = sub('beta','pval',betas)
  pvals.p = apply(pvals,2, function(y) paste0('P= ',formatC(y,format='g',digit=2)))
  if(nrow(res)==1) pvals.p = data.frame(t(pvals.p)) 
  colnames(pvals.p) = paste0('P',colnames(pvals.p))
  res = data.frame(res,ORs,pvals,pvals.p,stringsAsFactors=F)
}


# i bet rnd2 means round
rnd2 <- function(x) {
  round(x, 2)
}


# define wrapper function here
create_stratified_or_table_yi <- function(dat, env, snps, covs, mod){
  
  s = snps
  e = env
  
  res.pool.un = res.pool.g = res.pool.e = NULL
  Ncaco = data.frame(snps=snps,matrix(NA,ncol=6*2))
  colnames(Ncaco) = c('snps',paste0('Co',1:2),paste0('Ca',1:2),
                      paste0('Co1',1:2),paste0('Ca1',1:2),
                      paste0('Co2',1:2),paste0('Ca2',1:2))
  rownames(Ncaco) = snps
  
  # prep data 
  tmp1 = dat[, c('outcome', env, covs)]
  tmp2 = dat[, grepl(snps, names(dat))] # if you're incorporating Ref/Alt there's no issue with non-unique SNPs.. 
  
  # .... new part, try to recode depending on direction of association
  tmp_model <- paste0("outcome ~ ", paste0(snps, "_dose"), "+", env, "+", paste0(covs, collapse = "+"))
  tmp_coeff <- coef(lm(tmp_model, data = dat))[2]
  if(tmp_coeff > 0) {
    names(tmp2) <- c("dosage", "p0", "p1", "p2")
  } else {
    names(tmp2) <- c("dosage", "p2", "p1", "p0")
  }
  # names(tmp2) <- c("dosage", "p0", "p1", "p2")
  Data = cbind(tmp1, tmp2) %>% 
    na.omit(.[, c(env,'outcome', covs, 'p1','p2','dosage')])
  
  #======================================================================#
  #---- Calculate counts for each cell ---------------- #
  Ncaco[s,c(paste0('Co',1:2),paste0('Ca',1:2))] = t(table(Data[,c('outcome',e)]))
  Ncaco[s,c(paste0('Co1',1:2),paste0('Ca1',1:2))] = c(t(tapply(Data$p1,Data[,c('outcome',e)],sum,na.rm=T)))
  Ncaco[s,c(paste0('Co2',1:2),paste0('Ca2',1:2))] = c(t(tapply(Data$p2,Data[,c('outcome',e)],sum,na.rm=T)))
  
  #======================================================================#
  #-- Fit unrestricted model -------- #
  Data[,env] = factor(Data[,env]) # E should be a factor
  tmp = as.vector(Model.all.new(Data,mod,env)) # let's try fitting this model
  res.pool.un = rbind(res.pool.un,data.frame(s,env,t(tmp$GE)))
  res.pool.e = rbind(res.pool.e,data.frame(s,env,t(tmp$E)))
  res.pool.g = rbind(res.pool.g,data.frame(s,env,t(tmp$G)))
  
  ## organize results ##
  elvl=c(0:1) ; glvl = c(0,1,2)
  
  colnames(res.pool.un) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
                            'beta0p1','beta0p2','se0p1','se0p2',
                            paste0('beta',elvl[-1],'p1'),paste0('se',elvl[-1],'p1'),
                            paste0('beta',elvl[-1],'p2'),paste0('se',elvl[-1],'p2'))
  res.pool.un = format.res(res.pool.un)
  
  ##== stratified by G results ##
  colnames(res.pool.g) = c('snp','env',paste0(rep(c('beta0','se0','beta1','se1','beta2','se2'),each=1),
                                              rep(elvl[-1],6)))
  res.pool.g = format.res(res.pool.g)
  
  ##== stratified by E results ##
  colnames(res.pool.e) = c('snp','env',paste0('beta1',elvl),paste0('se1',elvl),
                           paste0('beta2',elvl),paste0('se2',elvl))
  res.pool.e = format.res(res.pool.e)
  
  ##== Put into table 
  OR.tab = ORtab(s,elvl=elvl,glvl=glvl,res=res.pool.un)
  pval.tab = ptab(s,elvl=elvl,glvl=glvl,res=res.pool.un)
  
  ORg.tab = matrix(as.character(unlist(list(res.pool.g[,paste0('OR0',elvl[-1])],
                                            res.pool.g[,paste0('OR1',elvl[-1])],
                                            res.pool.g[,paste0('OR2',elvl[-1])]) ) ),ncol=3)
  rownames(ORg.tab) <- 'OR2'
  pg.tab = matrix(as.character(unlist(list(res.pool.g[,paste0('Ppval0',elvl[-1])],
                                           res.pool.g[,paste0('Ppval1',elvl[-1])],
                                           res.pool.g[,paste0('Ppval2',elvl[-1])]))),ncol=3)
  
  colnames(ORg.tab) = colnames(pg.tab) = c(paste0('p',glvl))
  ORe.tab = matrix(as.character(unlist(c(res.pool.e[,paste0('OR1',elvl)],
                                         res.pool.e[,paste0('OR2',elvl)]))),ncol=2)
  pe.tab  = matrix(as.character(unlist(c(res.pool.e[,paste0('Ppval1',elvl)],
                                         res.pool.e[,paste0('Ppval2',elvl)]))),ncol=2,byrow=T)
  
  #== calculate counts for G=0 and put counts into format ca/co
  for(i in 1:2){
    eval(parse(text=paste0('Ncaco$caco0',i,"=paste0(round(Ncaco$Ca",i,'-Ncaco$Ca1',i,'-Ncaco$Ca2',i,
                           ",1),'/',round(Ncaco$Co",i,'-Ncaco$Co1',i,'-Ncaco$Co2',i,',1))')))
    eval(parse(text=paste0('Ncaco$caco1',i,'=paste0(round(Ncaco$Ca1',i,",1),'/',round(Ncaco$Co1",i,",1))")))
    eval(parse(text=paste0('Ncaco$caco2',i,'=paste0(round(Ncaco$Ca2',i,",1),'/',round(Ncaco$Co2",i,",1))")))
  }
  
  #=== write into table
  est = cbind(rbind(OR.tab[1,],pval.tab[1,],OR.tab[2,],pval.tab[2,],
                    ORg.tab[1,],pg.tab[1,]),
              rbind(ORe.tab[1,],pe.tab[1,],ORe.tab[2,],pe.tab[2,],rep(NA,2),rep(NA,2)),
              N0 = c(Ncaco[s,'caco01'],NA,Ncaco[s,'caco02'],rep(NA,3)),
              N1 = c(Ncaco[s,'caco11'],NA,Ncaco[s,'caco12'],rep(NA,3)),
              N2 = c(Ncaco[s,'caco21'],NA,Ncaco[s,'caco22'],rep(NA,3)))
  
  # class(est)
  # est
  snp_rename <- (gsub('\\.', '_', snps) %>%  gsub('X', 'chr', .))
  saveRDS(est, file = paste0("/media/work/tmp/posthoc/", env, "/stratified_or_dataframe_", env, "_", snp_rename, "_", paste0(sort(covs), collapse = "_"), ".rds"), version = 2)
  
}
  

# covs = c('age_ref_imp','study_gxe','pc1','pc2','pc3')
# mod = 'age_ref_imp+study_gxe+pc1+pc2+pc3+study_gxe:pc1+study_gxe:pc2+study_gxe:pc3+study_gxe:age_ref_imp+study_gxe:p1 + study_gxe:p2+ study_gxe:hrt_ref_pm2'
# env ='hrt_ref_pm2' ## E
# snps='X12.111401316.G.T' ##  SNP 
# s='X12.111401316.G.T' ##  SNP 
# 
# 
# 
# test <- create_stratified_or_table_yi(dat = dat, env = 'hrt_ref_pm2', snps='X12.111401316.G.T', covs = c('age_ref_imp','study_gxe','pc1','pc2','pc3'), mod = 'age_ref_imp+study_gxe+pc1+pc2+pc3+study_gxe:pc1+study_gxe:pc2+study_gxe:pc3+study_gxe:age_ref_imp+study_gxe:p1 + study_gxe:p2+ study_gxe:hrt_ref_pm2' )
# 
# test


# ---- Let's try applying function over a vector of SNPs ---- #
geno <- readRDS(paste0("/media/work/tmp/posthoc/gwis_sig_results_output_", env, ".rds"))
dat <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", env, "_basic_covars_glm.rds")) %>% 
  inner_join(geno, 'vcfid') %>% 
  rename_at(vars(c("PC1", "PC2", "PC3")), tolower) # for neatness

# flip E coding to risk increase
if(env %in% c("hrt_ref_pm2", "ep_ref_pm_gxe", "eo_ref_pm_gxe")) {
  # dat <- mutate(dat, !!env := abs(!!(sym(env)) - 1))
  dat[, env] <- abs(dat[, env] - 1)
}

# create vector of snp variables to apply wrapper over
snp_list <- names(dat)[grepl("X\\d{1,2}", names(dat))] # anything starting with a X# or X##
snp_list <- snp_list[grepl("?dose", snp_list)] # get dose only (shortcut because regex no good)
snp_list <- gsub("_dose", "", snp_list) # take out _dose
snp_list

# run wrapper function (it saves RDS objects to /media/work/tmp/posthoc/ENV)
# mod... needs to be specified in function call... 
# 
# sapply(snp_list, function(x) create_stratified_or_table_yi(dat = dat, env = env, snps=x, covs = covs, mod = 'age_ref_imp+study_gxe+pc1+pc2+pc3+study_gxe:pc1+study_gxe:pc2+study_gxe:pc3+study_gxe:age_ref_imp+study_gxe:p1 + study_gxe:p2+ study_gxe:hrt_ref_pm2' ))
sapply(snp_list[1], function(x) create_stratified_or_table_yi(dat = dat, env = env, snps=x, covs = covs, mod = mod))

