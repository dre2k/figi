#=============================================================================#
# FIGI GWIS posthoc analysis
# Updated 05/01/2020
#
# output model summaries (data.frame, HTML/stargazer) of SNP level GxE:
# - stratified by Sex, study-design, and other variables of interest
#   - output is an html file that contains stargazer code for Rmarkdown
#   - output saved under /media/work/gwis/posthoc/exposure
#
# IN PROGRESS
# - adjusted by additional covariates as required/requested
# - stratified odds ratios by SNP x E (Yi's code)
#
#
#
# need to specify model, I might limit this to basic covariates only
# extend later if requestd.......
#=============================================================================#
library(tidyverse)
library(data.table)
library(stargazer)
library(lmtest)
library(figifs)
rm(list = ls())

# exposure <- 'asp_ref'
# covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
# hrc_version <- 'v2.4'
# snp <- '1:15629040:T:G'
# snp <- paste0("chr", gsub("\\:", "\\_", snp))

# make sure order of arguments match
# args <- commandArgs(trailingOnly=T)
# exposure <- args[1] # ex: asp_ref
# hrc_version <- args[2] # ex: v2.4
# snp <- args[3] # ex: 2:27730940:T:C	
# snp <- paste0("chr", gsub("\\:", "\\_", snp))
# stat <- args[4]
# covariates <- sort(c(args[5:length(args)]))

# output_dir <- paste0("/media/work/gwis/posthoc/", exposure, "/")

# source("/home/rak/Dropbox/FIGI/FIGI_code/results/functions_ver2.R")
# dir.create(paste0("/media/work/gwis/posthoc/", exposure), showWarnings = FALSE)


#---------------------------#
# this might need to be a different script for convenience
#----------------------------#
# Yi code (stratified odds ratios)
# source("~/Dropbox/FIGI/FIGI_code/results/functions_yi.R")



# ----- dataset + genotypes ----- #
# vcfid_list <- readRDS(paste0("~/data/GxEScanR_PhenotypeFiles/FIGI_v2.3_gxeset_", exposure, "_basic_covars_glm.rds"))[,'vcfid']
vcfid_list <- readRDS(paste0("~/data/results/input/FIGI_", hrc_version , "_gxeset_", exposure, "_basic_covars_glm.rds"))[,'vcfid']
geno <- readRDS(paste0("/media/work/gwis/posthoc/gwis_sig_results_output_", exposure, ".rds"))

# prettify SNP names
tmp <- names(geno)
tmp <- gsub("\\.", "\\_", tmp)
tmp <- gsub("X", "chr", tmp)
# tmp <- gsub("\\_dose", "", tmp)
names(geno) <- tmp

# read in data, join, redefine some variables
# data already has study_design information
figi <- readRDS(paste0("~/data/results/input/FIGI_", hrc_version, "_gxeset_analysis_data_glm.rds")) %>%
  inner_join(geno, 'vcfid') %>%
  filter(vcfid %in% vcfid_list) %>%
  mutate(cancer_site_sum2 = factor(cancer_site_sum2, levels = c("proximal", "distal", "rectal")),
         vigmodlns = as.numeric(vigmodlns))

# ---- > this is only necessary if creating stratified odds ratios
# check direction of SNP, switch if necessary
snp_model <- lm(as.formula(paste0("outcome ~ ", snp, "_dose")), data = figi)

# flip dosages
if(snp_model$coefficients[2] < 0) {
  figi[, paste0(snp, '_dose')] <- abs(figi[, paste0(snp, '_dose')] - 2)
}

# flip genotype probabilities
if(snp_model$coefficients[2] < 0) {
  pp <- figi[,paste0(snp, "_p2")]
  figi[,paste0(snp, "_p2")] <- figi[, paste0(snp, "_p0")]
  figi[,paste0(snp, "_p0")] <- pp
}


tmp1 = figi[, c('outcome', exposure, covariates)]
tmp2 = figi[, grepl(snp, names(figi))]
# for linear regression and checking direction of SNP risk
tmp_model <- paste0("outcome ~ ", paste0(snp, "_dose"), "+", exposure, "+", paste0(covariates, collapse = "+"))
tmp_coeff <- coef(lm(tmp_model, data = figi))[2]
if(tmp_coeff > 0) {
  names(tmp2) <- c("dosage", "p0", "p1", "p2")
} else {
  names(tmp2) <- c("dosage", "p2", "p1", "p0")
}
Data = cbind(tmp1, tmp2) %>% 
  na.omit(.[, c(env,'outcome', covs, 'p1','p2','dosage')])










# covariates <- readRDS(paste0(output_dir, 'posthoc_helper_', exposure, '_basic_covariates.rds'))
# mod = readRDS(paste0(output_dir, 'posthoc_helper_', exposure, '_stratified_or_model.rds'))

mod = 'age_ref_imp+sex+pc1+pc2+pc3+study_gxe'




res.pool.un = res.pool.g = res.pool.e = NULL
Ncaco = data.frame(snps=snp,matrix(NA,ncol=6*2))
colnames(Ncaco) = c('snps',paste0('Co',1:2),paste0('Ca',1:2),
                    paste0('Co1',1:2),paste0('Ca1',1:2),
                    paste0('Co2',1:2),paste0('Ca2',1:2))
rownames(Ncaco) = snp



# Data = data[,c(1:12,grep(s,colnames(data)))]
# colnames(Data)[(ncol(Data)-2):ncol(Data)]=c('dosage','p1','p2')
# Data = na.omit(Data[,c(env,'outcome',covs,'p1','p2','dosage')])
#======================================================================
#---- Calculate counts for each cell ----------------
Ncaco[snp,c(paste0('Co',1:2),paste0('Ca',1:2))] = t(table(Data[,c('outcome',exposure)]))
Ncaco[snp,c(paste0('Co1',1:2),paste0('Ca1',1:2))] = c(t(tapply(Data$p1,Data[,c('outcome',exposure)],sum,na.rm=T)))
Ncaco[snp,c(paste0('Co2',1:2),paste0('Ca2',1:2))] = c(t(tapply(Data$p2,Data[,c('outcome',exposure)],sum,na.rm=T)))
#======================================================================
#-- Fit unrestricted model --------
Data[,exposure] = factor(Data[,exposure])
tmp = as.vector(Model.all.new(Data,mod,exposure))
res.pool.un = rbind(res.pool.un,data.frame(snp,exposure,t(tmp$GE)))
res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
res.pool.g = rbind(res.pool.g,data.frame(snp,exposure,t(tmp$G)))

## organize results ######
elvl=c(0:1) ; glvl = c(0,1,2)

# elvl=c(0:3) ; glvl = c(0,1,2)



colnames(res.pool.un) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
                          'beta0p1','beta0p2','se0p1','se0p2',
                          paste0('beta',elvl[-1],'p1'),paste0('se',elvl[-1],'p1'),
                          paste0('beta',elvl[-1],'p2'),paste0('se',elvl[-1],'p2'))
res.pool.un = format_res(res.pool.un)


##== stratified by G results #######
colnames(res.pool.g) = c('snp','env',paste0(rep(c('beta0','se0','beta1','se1','beta2','se2'),each=1),
                                            rep(elvl[-1],6)))
res.pool.g = format_res(res.pool.g)

##== stratified by E results #######
colnames(res.pool.e) = c('snp','env',paste0('beta1',elvl),paste0('se1',elvl),
                         paste0('beta2',elvl),paste0('se2',elvl))
res.pool.e = format_res(res.pool.e)

##== Put into table 
OR.tab = ORtab(snp,elvl=elvl,glvl=glvl,res=res.pool.un)
pval.tab = ptab(snp,elvl=elvl,glvl=glvl,res=res.pool.un)

ORg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('OR0',elvl[-1])]),
                                       as.character(res.pool.g[,paste0('OR1',elvl[-1])]),
                                       as.character(res.pool.g[,paste0('OR2',elvl[-1])])))), ncol=3)
pg.tab = matrix(as.character(unlist(c(as.character(res.pool.g[,paste0('Ppval0',elvl[-1])]),
                                      as.character(res.pool.g[,paste0('Ppval1',elvl[-1])]),
                                      as.character(res.pool.g[,paste0('Ppval2',elvl[-1])])))),ncol=3)

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
            N0 = c(Ncaco[snp,'caco01'],NA,Ncaco[snp,'caco02'],rep(NA,3)),
            N1 = c(Ncaco[snp,'caco11'],NA,Ncaco[snp,'caco12'],rep(NA,3)),
            N2 = c(Ncaco[snp,'caco21'],NA,Ncaco[snp,'caco22'],rep(NA,3)))




# and then generate that other column


tmp = as.vector(Model.all.new.dosage(Data,mod,exposure))
res.pool.un = rbind(res.pool.un,data.frame(snp,exposure,t(tmp$GE)))
res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
res.pool.g = rbind(res.pool.g,data.frame(snp,exposure,t(tmp$G)))



elvl=c(0:1) ; glvl = c(0,1)
test = NULL
test = rbind(test,data.frame(snp,exposure,t(tmp$GE)))

colnames(test) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
                   'beta0p1','beta0p2','se0p1','se0p2')
test = format_res(test)

res.pool.e = NULL
res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(huh$E)))
colnames(res.pool.e) = c('snp',
                         'env',
                         paste0('beta1',elvl),paste0('se1',elvl))
test <- format_res(res.pool.e)


ORe.tab = matrix(as.character(unlist(c(test[,paste0('OR1',elvl)]))),ncol=2)
pe.tab  = matrix(as.character(unlist(c(test[,paste0('Ppval1',elvl)]))),ncol=2,byrow=T)
est2 <- rbind(ORe.tab[1,1],pe.tab[1,1],ORe.tab[1,2],pe.tab[1,2], NA, NA)

final_out <- est %>% 
  dplyr::select(-c(4,5)) %>% 
  add_column(est2, .after = 3)
final_out <- add_column(est, est2, .after = 4)
