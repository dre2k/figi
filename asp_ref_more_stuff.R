# write function to organize the results for stratified odds ratios for manuscript (NSAIDs)

snp = "chr6_12577203_T_C"
covariates


stratified_or_organize <- function(exposure, snp, covariates) {
  tmp <- readRDS(glue("/media/work/gwis_test/{exposure}/output/posthoc/stratified_oddsratio_{exposure}_v2.4_{snp}_{glue_collapse(sort(covariates), sep = '_')}.rds")) %>% 
    setNames(make.unique(names(.))) %>% 
    separate(N0, into = c("N0_case", "N0_control"), sep = "/") %>% 
    separate(N1, into = c("N1_case", "N1_control"), sep = "/") %>% 
    separate(N2, into = c("N2_case", "N2_control"), sep = "/") %>% 
    mutate(across(.cols = starts_with("N"), ~ round(as.numeric(.x), 0) ), 
           across(.cols = where(is.numeric), ~ format(.x,big.mark=",",scientific=FALSE)))
  
  out <- matrix(nrow = 3, ncol = 9)
  out[1,] <- c(tmp$N0_case[1], tmp$N0_control[1], "1.0 (Ref)",
               tmp$N1_case[1], tmp$N1_control[1], "1.0 (Ref)",
               tmp$N2_case[1], tmp$N2_control[1], "1.0 (Ref)")
  out[2,] <- c(tmp$N0_case[3], tmp$N0_control[3], tmp[5,1],
               tmp$N1_case[3], tmp$N1_control[3], tmp[5,2],
               tmp$N2_case[3], tmp$N2_control[3], tmp[5,3])
  out[3,1] <- exposure
  out[3,2] <- snp
  out[3,3] <- glue_collapse(covariates, sep = "_")
  
  write.csv(out, "~/Dropbox/tmp.csv", quote = T, row.names = F, append = T)
}


stratified_or_organize(exposure, "chr6_12577203_T_C", covariates)
covariates_mult = c(covariates, "smk_ever", "bmi", "alcoholc", "redmeatqc2")

stratified_or_organize(exposure, "chr6_12577203_T_C", covariates_mult)



stratified_or_organize(exposure, "chr5_40252294_C_T", covariates)
covariates_mult = c(covariates, "smk_ever", "bmi", "alcoholc", "redmeatqc2")

stratified_or_organize(exposure, "chr5_40252294_C_T", covariates_mult)










