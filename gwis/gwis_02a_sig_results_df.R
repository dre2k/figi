# --------------------------------------------------------------------------- #
# gwis results for G, GxE, E|G (case, control, combined), 2DF, 3DF
# assemble significant hits into a data.frame for presentation in Rmarkdown
# for some methods, there are many SNPs that meet sig threshold (e.g. 2df or suggestive gxe hits)
# clump these results, and present the top SNPs in each clump for the table of results
#
# with and without exclusion of GWAS, D|G, and E|G (controls) hits for 2DF/3DF
# 
# note that some variables are defined globally in parent script (gwis_03...) 
# --------------------------------------------------------------------------- #

# calculate chisq statistic p values
# output SNP, Chr, BP, Ref/Alt, N, Cases, and p-value
plot_sig_wrapper <- function(dat, statistic, df, sig_level, hrc_version = hrc_version, exclude = "") {
  
  # exclude gwas loci or SNPs with strong marginal effects
  if(exclude == '_no_gwas') {
    tmp <- dat %>%
      dplyr::filter(!SNP2 %in% exclude_gwas)
  } else if(exclude == '_no_gwas_no_marginal') {
    tmp <- dat %>%
      dplyr::filter(!SNP2 %in% exclude_gwas,
                    !SNP %in% exclude_g)
  } else if(exclude == '_no_gwas_no_marginal_no_ge') {
    tmp <- dat %>%
      dplyr::filter(!SNP2 %in% exclude_gwas,
                    !SNP %in% exclude_g,
                    !SNP %in% exclude_eg_control)
  } else if(exclude == '_no_gwas_no_ge') {
    tmp <- dat %>%
      dplyr::filter(!SNP2 %in% exclude_gwas,
                    !SNP %in% exclude_eg_control)
  } else {
    tmp <- dat
  }
  
  # calculate pvalue, filter p value by significance level, output specific columns
  tmp <- tmp %>%
    mutate(Pval = calculate_pval((!!sym(statistic)), df)) %>%
    dplyr::filter(Pval < sig_level) %>%
    mutate(Pval = formatC(Pval, format = 'e', digits = 2)) %>%
    dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval)
  
  saveRDS(tmp, file = paste0(output_dir, "significant_results_dataframe_", statistic, "_", exposure, exclude, ".rds"), version = 2)
}

# ------- no exclusions ------- #
# Marginal G Results
plot_sig_wrapper(gxe, 'chiSqG', df = 1, sig_level = 5e-8)
# GxE results - p<5e-6
plot_sig_wrapper(gxe, 'chiSqGxE', df = 1, sig_level = 5e-6)
# 2DF results
plot_sig_wrapper(gxe, 'chiSq2df', df = 2, sig_level = 5e-8)
# 3DF results
plot_sig_wrapper(gxe, 'chiSq3df', df = 3, sig_level = 5e-8)
# GE
plot_sig_wrapper(gxe, 'chiSqGE', df = 1, sig_level = 5e-8)
# GE Among Controls
plot_sig_wrapper(gxe, 'chiSqControl', df = 1, sig_level = 5e-8)
# GE Case Only
plot_sig_wrapper(gxe, 'chiSqCase', df = 1, sig_level = 5e-8)



# -------- exclude GWAS, D|G, and E|G (controls) hits from 2DF/3DF results ------- #
# wrapper function
# reformat table to include inidividual components of joint tests
plot_sig_wrapper <- function(dat, statistic, df, sig_level, hrc_version = hrc_version, exclude = ""){
  
  # exclude gwas loci or SNPs with strong marginal effects
  if(exclude == '_no_gwas') {
    tmp <- dat %>%
      dplyr::filter(!SNP2 %in% exclude_gwas)
  } else if(exclude == '_no_gwas_no_marginal') {
    tmp <- dat %>%
      dplyr::filter(!SNP2 %in% exclude_gwas,
                    !SNP %in% exclude_g)
  } else if(exclude == '_no_gwas_no_marginal_no_ge') {
    tmp <- dat %>%
      dplyr::filter(!SNP2 %in% exclude_gwas,
                    !SNP %in% exclude_g,
                    !SNP %in% exclude_eg_control)
  } else if(exclude == '_no_gwas_no_ge') {
    tmp <- dat %>%
      dplyr::filter(!SNP2 %in% exclude_gwas,
                    !SNP %in% exclude_eg_control)
  } else {
    tmp <- dat
  }
  
  if(statistic == 'chiSq2df') {
    tmp <- tmp %>%
      mutate(Pval_2df = calculate_pval((!!sym(statistic)), df),
             Pval_dg = calculate_pval(chiSqG, df = 1), 
             Pval_gxe = calculate_pval(chiSqGxE, df = 1), 
             Pval_eg = calculate_pval(chiSqGE, df = 1)) %>%
      dplyr::filter(Pval_2df < sig_level) %>%
      mutate(Pval_2df = formatC(Pval_2df, format = 'e', digits = 2),
             Pval_dg = formatC(Pval_dg, format = 'e', digits = 2),
             Pval_gxe = formatC(Pval_gxe, format = 'e', digits = 2)) %>%
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval_dg, Pval_gxe, Pval_2df)
  } else if(statistic == 'chiSq3df') {
    tmp <- tmp %>%
      mutate(Pval_3df = calculate_pval((!!sym(statistic)), df),
             Pval_dg = calculate_pval(chiSqG, df = 1), 
             Pval_gxe = calculate_pval(chiSqGxE, df = 1), 
             Pval_eg = calculate_pval(chiSqGE, df = 1)) %>%
      dplyr::filter(Pval_3df < sig_level) %>%
      mutate(Pval_3df = formatC(Pval_3df, format = 'e', digits = 2), 
             Pval_dg = formatC(Pval_dg, format = 'e', digits = 2),
             Pval_gxe = formatC(Pval_gxe, format = 'e', digits = 2),
             Pval_eg = formatC(Pval_eg, format = 'e', digits = 2)) %>%
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval_dg, Pval_gxe, Pval_eg, Pval_3df)
  }
  
  saveRDS(tmp, file = paste0(output_dir, "significant_results_dataframe_", statistic, "_", exposure, exclude, ".rds"), version = 2)
}




# ------- exclude known loci ------- #
# 2DF results
plot_sig_wrapper(gxe, 'chiSq2df', df = 2, sig_level = 5e-8, exclude = '_no_gwas')
# 3DF results
plot_sig_wrapper(gxe, 'chiSq3df', df = 3, sig_level = 5e-8, exclude = '_no_gwas')
# 2DF results
plot_sig_wrapper(gxe, 'chiSq2df', df = 2, sig_level = 5e-8, exclude = '_no_gwas_no_marginal')
# 3DF results
plot_sig_wrapper(gxe, 'chiSq3df', df = 3, sig_level = 5e-8, exclude = '_no_gwas_no_marginal')
# 3DF results
plot_sig_wrapper(gxe, 'chiSq3df', df = 3, sig_level = 5e-8, exclude = '_no_gwas_no_marginal_no_ge')
# 3DF results
plot_sig_wrapper(gxe, 'chiSq3df', df = 3, sig_level = 5e-8, exclude = '_no_gwas_no_ge')



# ------- clumped outputs -------- #
plot_sig_wrapper_clump <- function(dat, statistic, df, sig_level, hrc_version = hrc_version, exclude = "_ld_clump_no_gwas") {

  tmp <- dat %>%
    dplyr::filter(!SNP2 %in% exclude_gwas)

  if(statistic == 'chiSq2df') {
    tmp <- tmp %>%
      mutate(Pval_2df = calculate_pval((!!sym(statistic)), df),
             Pval_dg = calculate_pval(chiSqG, df = 1),
             Pval_gxe = calculate_pval(chiSqGxE, df = 1),
             Pval_eg = calculate_pval(chiSqGE, df = 1)) %>%
      dplyr::filter(Pval_2df < sig_level) %>%
      mutate(Pval_2df = formatC(Pval_2df, format = 'e', digits = 2),
             Pval_dg = formatC(Pval_dg, format = 'e', digits = 2),
             Pval_gxe = formatC(Pval_gxe, format = 'e', digits = 2)) %>%
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval_dg, Pval_gxe, Pval_2df)
  } else if(statistic == 'chiSq3df') {
    tmp <- tmp %>%
      mutate(Pval_3df = calculate_pval((!!sym(statistic)), df),
             Pval_dg = calculate_pval(chiSqG, df = 1),
             Pval_gxe = calculate_pval(chiSqGxE, df = 1),
             Pval_eg = calculate_pval(chiSqGE, df = 1)) %>%
      dplyr::filter(Pval_3df < sig_level) %>%
      mutate(Pval_3df = formatC(Pval_3df, format = 'e', digits = 2),
             Pval_dg = formatC(Pval_dg, format = 'e', digits = 2),
             Pval_gxe = formatC(Pval_gxe, format = 'e', digits = 2),
             Pval_eg = formatC(Pval_eg, format = 'e', digits = 2)) %>%
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval_dg, Pval_gxe, Pval_eg, Pval_3df)
  } else {
    tmp <- tmp %>%
      mutate(Pval = calculate_pval((!!sym(statistic)), df)) %>%
      dplyr::filter(Pval < sig_level) %>%
      mutate(Pval = formatC(Pval, format = 'e', digits = 2)) %>%
      dplyr::select(SNP, Chromosome, Location, Reference, Alternate, Subjects, Cases, Pval)
  }

  saveRDS(tmp, file = paste0(output_dir, "significant_results_dataframe_", statistic, "_", exposure, exclude, ".rds"), version = 2)
}


# GxE

# NEED TO FILTER OUT CLUMPED SNPS USING CLUMPED OUTPUT!

# ------- no exclusions ------- #
# GxE results - p<5e-6
plot_sig_wrapper_clump(gxe, 'chiSqGxE', df = 1, sig_level = 5e-6)
# 2df
plot_sig_wrapper_clump(gxe, 'chiSq2df', df = 2, sig_level = 5e-8)
# 3df
plot_sig_wrapper_clump(gxe, 'chiSq3df', df = 3, sig_level = 5e-8)
