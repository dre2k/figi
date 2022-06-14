# generating additional analyses for aspirin paper for god's sake



data_epi <- input_data
exposure 
snp <- '5:40252294:C:T'
covariates

method = 'chiSqGxE'
flip_allele = F
path = '/media/work/gwis_test/asp_ref/'
path = '/media/work/gwis_test/aspirin/'



fit_gxe_stratified <- function(data_epi,
                               exposure,
                               snp,
                               covariates,
                               strata = c('sex', 'study_design', 'cancer_site_sum2', 'bmic3'),
                               method = c('chiSqGxE', 'chiSqCase', 'chiSq2df', 'chiSq3df', 'chiSqGE'),
                               flip_allele = F,
                               path) {
  
  # make path the exposure directory e.g. /media/work/gwis_test/exposure
  wdir = glue("{path}/output/posthoc")
  
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
  
  data_dose <-
    qread(glue("{wdir}/dosage_chr{snp_info[1]}_{snp_info[2]}.qs"))
  data <- inner_join(data_epi, data_dose, 'vcfid')
  snp_new <- snpfix

  # let's function and just output the data.frame to make things easier
  

  s1a <- which(data[, 'sex'] == 0 & (data[, 'cancer_site_sum2'] == "proximal" | data[, 'outcome'] == 0))
  m1a <- fit_gxe(ds = data[s1a, ], exposure = exposure, snp = 'chr5_40252294_C_T_dose', covariates = c("age_ref_imp", "pc1", "pc2", "pc3", "study_gxe"), method = 'chiSqGxE') 
  
  
  s2a <- which(data[, 'sex'] == 0 & (data[, 'cancer_site_sum2'] == "distal" | data[, 'outcome'] == 0))
  m2a <- fit_gxe(ds = data[s2a, ], exposure = exposure, snp = 'chr5_40252294_C_T_dose', covariates = c("age_ref_imp", "pc1", "pc2", "pc3", "study_gxe"), method = 'chiSqGxE') 
  
  
  s3a <- which(data[, 'sex'] == 0 & (data[, 'cancer_site_sum2'] == "rectal" | data[, 'outcome'] == 0))
  m3a <- fit_gxe(ds = data[s3a, ], exposure = exposure, snp = 'chr5_40252294_C_T_dose', covariates = c("age_ref_imp", "pc1", "pc2", "pc3", "study_gxe"), method = 'chiSqGxE') 
  
  s1b <- which(data[, 'sex'] == 1 & (data[, 'cancer_site_sum2'] == "proximal" | data[, 'outcome'] == 0))
  m1b <- fit_gxe(ds = data[s1b, ], exposure = exposure, snp = 'chr5_40252294_C_T_dose', covariates = c("age_ref_imp", "pc1", "pc2", "pc3", "study_gxe"), method = 'chiSqGxE') 
  
  
  s2b <- which(data[, 'sex'] == 1 & (data[, 'cancer_site_sum2'] == "distal" | data[, 'outcome'] == 0))
  m2b <- fit_gxe(ds = data[s2b, ], exposure = exposure, snp = 'chr5_40252294_C_T_dose', covariates = c("age_ref_imp", "pc1", "pc2", "pc3", "study_gxe"), method = 'chiSqGxE') 
  
  
  s3b <- which(data[, 'sex'] == 1 & (data[, 'cancer_site_sum2'] == "rectal" | data[, 'outcome'] == 0))
  m3b <- fit_gxe(ds = data[s3b, ], exposure = exposure, snp = 'chr5_40252294_C_T_dose', covariates = c("age_ref_imp", "pc1", "pc2", "pc3", "study_gxe"), method = 'chiSqGxE') 
  
  
  

out <- list(m1a, m2a, m3a, m1b, m2b, m3b)

out_df <- do.call(bind_rows, map(out, ~ tidy(.x[[1]])))
write.csv(out_df, file = "~/Dropbox/aspirin_chr5_sex_tumorsite_stratified.csv", quote = F, row.names = F)

