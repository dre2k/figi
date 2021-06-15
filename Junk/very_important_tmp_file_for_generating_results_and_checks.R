
# =========================================================================== #
# posthoc function wrappers
# =========================================================================== #

# pvalue table, gxe models
posthoc_run_models <- function(x, method) {
  walk(x, ~ pval_summary(ds = figi, exposure, .x, covariates_list, method, output_dir = output_dir))
  walk(x, ~ fit_gxe_covars(ds = figi, exposure, .x, covariates_list, method, output_dir = output_dir))
  walk(x, ~ fit_stratified_or(ds = figi, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates_set1, mod = mod, dosage = F, output_dir = output_dir))
}


# reri and interaction plot
posthoc_run_reri_plot <- function(y) {
  # RERI
  walk(y, ~ reri_wrapper(x = figi, exposure = exposure, snp = .x, covariates = covariates, output_dir = output_dir))
  # interaction plot
  walk(y, ~ iplot_wrapper(x = figi, exposure = exposure, snp = .x, covariates = covariates, output_dir = output_dir)) 
}


# locuszoom and functional annotation plots
posthoc_create_plots <- function(x, statistic_to_plot) {
  
  locuszoom_helper <- function(snp) {
    tmp <- gsub("chr", "", snp) # have to turn snps into chr:bp format for locuszoom 
    snp_chr <- strsplit(tmp, "_")[[1]][1]
    snp_bp <- as.numeric(strsplit(tmp, '_')[[1]][2])
    paste(snp_chr, snp_bp, sep = ":")
  }
  
  # --- going to cat empty annotation to the denote-markers-file to avoid errors --- #
  snps_plot <- map_chr(x, locuszoom_helper)
  
  tmp <- read.table("/media/work/gwis/locuszoom/locuszoom_gecco_gwas_annotation.txt", header = T)
  
  locuszoom_marker_helper <- function(x) {
    tmp_xx <- c(paste0("chr", x), "", "white")
    tmp_out <- rbind(tmp, tmp_xx)
    write.table(tmp_out, file = paste0("/media/work/gwis/locuszoom/locuszoom_annotation_chr", x, "_", statistic_to_plot, ".txt"), quote = F, row.names = F, sep = '\t')
  }
  
  walk(snps_plot, locuszoom_marker_helper)
  
  walk(snps_plot, ~ system(paste("bash ~/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_02_locuszoom.sh", exposure, hrc_version, .x, statistic_to_plot)))
  walk(snps_plot, ~ system(paste("Rscript ~/Dropbox/FIGI/FIGI_code/results/posthoc/posthoc_03_functional_annotation.R", exposure, .x, statistic_to_plot)))
}










# 
# pvalue table, gxe models
posthoc_run_models <- function(x, method) {
  walk(x, ~ pval_summary(ds = figi, exposure, .x, covariates_list, method, output_dir = output_dir))
  walk(x, ~ fit_gxe_covars(ds = figi, exposure, .x, covariates_list, method, output_dir = output_dir))
  walk(x, ~ fit_stratified_or_q4(ds = figi, exposure = exposure, snp = .x, hrc_version = hrc_version, covariates = covariates_set1, mod = mod, dosage = F, output_dir = output_dir))
}





#' fit_stratified_or
#' 
#' create stratified odds ratios table (Yi Lin)
#'
#' @param ds dataset
#' @param exposure string; exposure variable
#' @param snp string; snp name (chr1_bp_ref_alt)
#' @param hrc_version string; e.g. v2.4
#' @param covariates vector; adjustment covarites
#' @param mod string; model of exposure + covariates only. can specify study_gxe interactions if needed (that's how Yi's original code modeled these interaction models)
#' @param dosage whether you want to report dosage or genotype probabilities
#' @param output_dir string output directory
#'
#' @return a matrix of odds ratios, p-values, and strata counts
#' @export
#'
#' @examples fit_stratified_or(figi, 'asp_ref', 'chr5_1234_A_C', 'v2.4', c("age_ref_imp", "sex"), "age_ref_imp+sex")
fit_stratified_or_tmp <- function(ds, exposure, snp, hrc_version, covariates, mod, dosage = F, output_dir) {
  
  
  # SNP information
  snp_info <- unlist(strsplit(snp, split = "_"))
  
  model_check <- glm(glue("outcome ~ {exposure}*{snp} + {glue_collapse(covariates, sep = '+')}"), family = 'binomial', data = ds)
  
  # ---- recode SNPs so that lower risk allele is reference (to match RERI output)
  if (model_check[[1]][snp] >= 0) {
    # flip dosages
    ds[, snp] <- abs(ds[, paste0(snp)] - 2)
    ds[, paste0(snp, '_dose')] <- abs(ds[, paste0(snp, '_dose')] - 2)
    # flip genotype probabilities
    pp <- ds[,paste0(snp, "_p2")]
    ds[,paste0(snp, "_p2")] <- ds[, paste0(snp, "_p0")]
    ds[,paste0(snp, "_p0")] <- pp
    # assign ref/alt allele
    ref_allele <- snp_info[4]
    alt_allele <- snp_info[3]
  } else {
    ref_allele <- snp_info[3]
    alt_allele <- snp_info[4]
  }
  
  
  # create data subset
  tmp1 = ds[, c('outcome', exposure, covariates)]
  # tmp2 = ds[, grepl(snp, names(ds))]
  tmp2 = ds[, grepl(paste(paste0(snp, c("_dose", "_p0", "_p1", "_p2")), collapse = "|"), names(ds))]
  names(tmp2) <- c("dosage", "p0", "p1", "p2")
  
  ds_tmp = cbind(tmp1, tmp2) %>% 
    na.omit(.[, c(exposure,'outcome', covariates, 'p1','p2','dosage')])
  
  res.pool.un = res.pool.g = res.pool.e = NULL
  Ncaco = data.frame(snps=snp,matrix(NA,ncol=6*2))
  colnames(Ncaco) = c('snps',paste0('Co',1:2),paste0('Ca',1:2),
                      paste0('Co1',1:2),paste0('Ca1',1:2),
                      paste0('Co2',1:2),paste0('Ca2',1:2))
  rownames(Ncaco) = snp
  
  
  #---- Calculate counts for each cell ----------------
  Ncaco[snp,c(paste0('Co',1:2),paste0('Ca',1:2))] = t(table(ds_tmp[,c('outcome',exposure)]))
  Ncaco[snp,c(paste0('Co1',1:2),paste0('Ca1',1:2))] = c(t(tapply(ds_tmp$p1, ds_tmp[,c('outcome',exposure)],sum,na.rm=T)))
  Ncaco[snp,c(paste0('Co2',1:2),paste0('Ca2',1:2))] = c(t(tapply(ds_tmp$p2, ds_tmp[,c('outcome',exposure)],sum,na.rm=T)))
  
  #-- Fit unrestricted model --------
  ds_tmp[,exposure] = factor(ds_tmp[,exposure])
  tmp = as.vector(Model.all.new(ds_tmp,mod,exposure))
  res.pool.un = rbind(res.pool.un,data.frame(snp,exposure,t(tmp$GE)))
  res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
  res.pool.g = rbind(res.pool.g,data.frame(snp,exposure,t(tmp$G)))
  
  ## organize results ######
  elvl=c(0:1) ; glvl = c(0,1,2)
  
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
  
  # slightly modify output in case anyone wants per allele effects rather than p1/p2
  if(dosage == T) {
    
    tmp = as.vector(Model.all.new.dosage(ds_tmp,mod,exposure))
    
    res.pool.un = res.pool.g = res.pool.e = NULL
    res.pool.un = rbind(res.pool.un,data.frame(snp,exposure,t(tmp$GE)))
    # res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
    # res.pool.g = rbind(res.pool.g,data.frame(snp,exposure,t(tmp$G)))
    
    elvl=c(0:1) ; glvl = c(0,1)
    tmp2 = NULL
    tmp2 = rbind(tmp2,data.frame(snp,exposure,t(tmp$GE)))
    
    colnames(tmp2) = c('snp','env',paste0('beta',elvl[-1],'p0'),paste0('se',elvl[-1],'p0'),
                       'beta0p1','beta0p2','se0p1','se0p2')
    tmp2 = format_res(tmp2)
    
    res.pool.e = rbind(res.pool.e,data.frame(snp,exposure,t(tmp$E)))
    colnames(res.pool.e) = c('snp',
                             'env',
                             paste0('beta1',elvl),paste0('se1',elvl))
    tmp2 <- format_res(res.pool.e)
    
    
    ORe.tab = matrix(as.character(unlist(c(tmp2[,paste0('OR1',elvl)]))),ncol=2)
    pe.tab  = matrix(as.character(unlist(c(tmp2[,paste0('Ppval1',elvl)]))),ncol=2,byrow=T)
    est2 <- rbind(ORe.tab[1,1],pe.tab[1,1],ORe.tab[1,2],pe.tab[1,2], NA, NA)
    
    final_out <- est %>% 
      dplyr::select(-c(4,5)) %>% 
      add_column(est2, .after = 3)
    
    colnames(final_out) <- c(paste0(ref_allele, ref_allele), 
                             paste0(ref_allele, alt_allele), 
                             paste0(alt_allele, alt_allele), 
                             paste0(alt_allele, " Allelic Dosage"), 
                             "N0", "N1", "N2")
    # this is only for FACTOR variables (might have to modify when you run Q4 variables)
    exposure_level <- levels(ds[,exposure])
    rownames(final_out) <- c(paste0(exposure,"=",exposure_level[1]), 
                             "p0",
                             paste0(exposure,"=",exposure_level[2]), 
                             "p1", 
                             exposure, 
                             "p")
    saveRDS(final_out, file = paste0(output_dir, "stratified_oddsratio_", snp, "_", exposure, "_dosage.rds"))
  } else {
    colnames(est) <- c(paste0(ref_allele, ref_allele), 
                       paste0(ref_allele, alt_allele), 
                       paste0(alt_allele, alt_allele), 
                       paste0(ref_allele, alt_allele), 
                       paste0(alt_allele, alt_allele),
                       "N0", "N1", "N2")
    exposure_level <- levels(ds[,exposure])
    rownames(est) <- c(paste0(exposure,"=",exposure_level[1]), 
                       "p0",
                       paste0(exposure,"=",exposure_level[2]), 
                       "p1", 
                       exposure, 
                       "p")
    
    saveRDS(est, file = paste0(output_dir, "stratified_oddsratio_", snp, "_", exposure, "_norecode.rds"))
  }
  
}





fit_gxe_covars_tmp <- function(ds, 
                           exposure, 
                           snp, 
                           covariates_list, 
                           method = c('chiSqGxE', 'two-step', 'chiSqCase', 'chiSq2df', 'chiSq3df'),
                           output_dir) {
  
  method <- match.arg(method)
  
  # we probably don't want to flip SNPs for each covariates set, so let's base direction on the first set (always the simpler model), in the interaction model. 
  model_check <- glm(glue("outcome ~ {exposure}*{snp} + {glue_collapse(covariates_list[[1]], sep = '+')}"), family = 'binomial', data = ds)
  
  snp_old <- snp
  snp_tmp <- strsplit(snp, split = "_")
  chr <- snp_tmp[[1]][1]
  bp <- snp_tmp[[1]][2]
  a1 <- snp_tmp[[1]][3]
  a2 <- snp_tmp[[1]][4]
  
  snp_new <- glue("{chr}_{bp}_{a2}_{a1}_flipped")
  ds[[snp_new]] <- abs(2-ds[, snp_old])
  ref_allele = a2

  
  # if (model_check[[1]][snp] <  0) {
  #   snp_new <- glue("{chr}_{bp}_{a2}_{a1}_flipped")
  #   ds[[snp_new]] <- abs(2-ds[, snp_old])
  #   ref_allele = a2
  # } else {
  #   snp_new <- snp
  #   ref_allele = a1
  # }
  
  # apply 'fit_gxe' over covariate_list
  out <- lapply(covariates_list, function(x) fit_gxe(ds, exposure, snp_new, covariates = x))
  
  # combine them to call stargazer
  list_of_glms <- lapply(out, function(x) x[[1]])
  list_of_samplesizes <- lapply(list_of_glms, function(x) paste0(c("Ca=", "Co="), rev(as.character(table(x$model$outcome))), collapse = ','))
  coefs <- lapply(list_of_glms, function(x) (exp(coef(x))))
  
  # need to calculate p values
  list_of_chisq <- lapply(out, function(x) x[[2]])
  col_label = paste0(paste0("Covariate Set ", seq(1, length(out))), " (", list_of_samplesizes, ")")
  
  if(method == "chiSqGxE" | method == "twostep" | method == "chiSqCase") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[1]], df = 1, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("Reference allele = ", ref_allele),
               paste0(col_label, ", LRtest GxE p = ", gxe_pvalues))
  } else if(method == "chiSq2df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq(x[[2]], df = 2, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("Reference allele = ", ref_allele),
               paste0(col_label, ", LRtest 2DF p = ", gxe_pvalues))
  } else if(method == "chiSq3df") {
    gxe_pvalues <- do.call(c, lapply(list_of_chisq, function(x) formatC(pchisq((x[[2]] + x[[3]]), df = 3, lower.tail = F), format = "e", digits = 5)))
    notes <- c("(PC and Study estimates omitted from table)", 
               paste0("Reference allele = ", ref_allele),
               paste0(col_label, ", LRtest 3DF p = ", gxe_pvalues))
  }
  
  # save output html from stargazer
  out_html <- stargazer_helper(list_of_glms,
                               title=paste0(gsub("\\_", "\\\\_", snp_new), " x ", gsub('\\_', '\\\\_', exposure)), 
                               column.labels=col_label,
                               coef=coefs, 
                               notes=notes, single.row = T)
  
  # write object to html file
  cat(paste(out_html, collapse = "\n"), "\n", 
      file = paste0(output_dir, "gxe_", method, "_", snp, "_", exposure, "_covariate_sets.html"), append = F)
}










# =========================================================================== #
# just checknig how to remove GWAS hits from expectation based hybrid plots
# =========================================================================== #

x <- readRDS("/media/work/gwis/posthoc/diab/expectation_hybrid/diab_snplist_twostep_chiSqG_bin1_output.rds")
names(x)

exclude_gwas <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt") %>% 
  mutate(snps = paste0("X", Chr, ".", Pos )) %>% 
  pull(snps)

xnames <- substr(names(x), 1, nchar(names(x)) - 4)
xnames

x_index <- !xnames %in% exclude_gwas



# remove columns based on vector 'exclude_gwas'

x_out <- x[, x_index ]





#------------ debug -----------#


x = gxe
exposure
covariates = covariates
simplem_step1_statistic = 'chiSqG'
output_dir
filename_suffix = ""
include_gwas = F



testing <- function(x, exposure, covariates, simplem_step1_statistic, output_dir, filename_suffix = "", include_gwas=T) {
  files_input <- mixedsort(list.files(paste0(output_dir, "expectation_hybrid"), pattern = paste0(paste0("twostep_", simplem_step1_statistic, "_bin"), "(?:.+)", "output.rds"), full.names = T))
  files_list <- map(files_input, ~ readRDS(.x))
  
  
  #----- delete if you f* it up -----#
  #
  exclude_gwas_snps <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt") %>% 
    mutate(snps = paste0("X", Chr, ".", Pos )) %>% 
    pull(snps)
  
  tmp_function <- function(zz) {
    zznames <- substr(names(zz), 1, nchar(names(zz)) - 4)
    zz_index <- !zznames %in% exclude_gwas_snps
    zz_out <- zz[, zz_index]
    return(zz_out)
  }
  
  if(include_gwas==F) {
    files_list <- map(files_list, ~ tmp_function(.x))
  }
  #
  #----- delete if you f* it up -----#
  
  
  number_of_snps <- map_int(files_list, ~ ncol(.x)) - 1 # -1 to remove vcfid column 
  number_of_tests <- map_int(files_list, ~ meff_r(dat = .x, PCA_cutoff = 0.995, fixLength = 150))

exclude_gwas_snps <- fread("~/data/Annotations/gwas_141_ld_annotation_july2020.txt") %>% 
  mutate(snps = paste0("X", Chr, ".", Pos )) %>% 
  pull(snps)

tmp_function <- function(zz) {
  zznames <- substr(names(zz), 1, nchar(names(zz)) - 4)
  zz_index <- !zznames %in% exclude_gwas_snps
  zz_out <- zz[, zz_index]
  return(zz_out)
}

if(include_gwas==F) {
  files_list <- map(files_input, ~ tmp_function(.x))
}





# =========================================================================== #
# locuszoom plot annotation spitting out error
# need to make sure at least one line in the annotation file is in region being plotted
# include a blank line centered at the SNP being plotted.. need to create 
# annotation files for each SNP (easiest way)
# =========================================================================== #

xx <- read.table("/media/work/gwis/locuszoom/locuszoom_gecco_gwas_annotation.txt", header = T)
xx_tmp <- c("chr12:12345", "", "white")


xx_out <- rbind(xx, xx_tmp)
write.table(xx_out, file = "/media/work/gwis/locuszoom/test3.txt", row.names = F, quote = F, sep = "\t")


1.26  - 1.01 - 1.52 + 1
2.62 - 1.51 - 1.49 + 1
