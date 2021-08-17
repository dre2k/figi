#=============================================================================#
# FIGI GxE asp_ref results
#=============================================================================#
library(tidyverse)
library(data.table)
library(qqman)
library(EasyStrata)
library(lmtest)
library(figifs)
library(forcats)
library(car)
library(grid)
library(gridExtra)
library(stargazer)
library(nnet)
library(glue)
library(ramwas)
library(flextable)
library(gtools)
library(interactionR)
library(epiR)
library(flextable)
library(jtools)
library(interactions)
library(msm)
library(qs)
rm(list = ls())

# input variables
exposure = 'asp_ref'
hrc_version = 'v2.4'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset)



#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13)


# ------- stratified pooled analysis ------- #
pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'sex', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_sex"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_glm(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, strata = 'study_design', filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_study_design"), output_dir = glue("{path}/output/posthoc/"))

pooled_analysis_multinom(input_data, exposure = exposure, hrc_version = hrc_version, covariates = covariates, filename_suffix = paste0(paste0(sort(covariates), collapse = '_'), "_stratified_cancer_site_sum2"), output_dir = glue("{path}/output/posthoc/"))





# main effects additional ----
output_dir_dropbox = paste0("~/Dropbox/Presentations/", exposure, "/")

tmp1 <- readRDS(paste0("/media/work/gwis/results/input/FIGI_", hrc_version, "_gxeset_", exposure, "_basic_covars_glm.rds")) %>% 
  pull(vcfid)
xx <- filter(input_data, vcfid %in% tmp1) %>% 
  group_by(study_gxe)


results_beta <- dplyr::do(xx, broom::tidy(glm(outcome ~ folate_dietqc2 + age_ref_imp + sex + energytot_imp, data = . , family = 'binomial'))) %>% 
  dplyr::filter(grepl("folate_dietqc2", term)) %>% 
  dplyr::arrange(study_gxe) %>% 
  inner_join(unique(xx[,c('study_gxe', 'study_design')]), 'study_gxe')

results_meta <- meta::metagen(estimate,
                              std.error,
                              data=results_beta,
                              studlab=paste(study_gxe),
                              comb.fixed = FALSE,
                              comb.random = TRUE,
                              method.tau = "SJ",
                              hakn = TRUE,
                              prediction=TRUE,
                              sm="OR", 
                              byvar=study_design)

fo <- find.outliers(results_meta)
fo
meta::forest(results_meta,
             layout = "JAMA",
             text.predict = "95% CI",
             col.predict = "black",
             # leftcols = c("studlab", "Control", "Case", "N", "effect", "ci", "w.random"),
             digits.addcols=0,
             study.results=T,
             prediction = F,
             col.random = 'red')

png(paste0(output_dir_dropbox, "meta_analysis_", "folate_dietqc2",  "_", "original_outliers_removed", ".png"), height = 17, width = 8.5, units = 'in', res = 150)                                                          
forest(fo)
dev.off()


# leave on out (influence analysis)
inf.analysis <- InfluenceAnalysis(x = results_meta,
                                  random = TRUE)

summary(inf.analysis)
plot(inf.analysis, "influence")
plot(inf.analysis, "baujat")
plot(inf.analysis, "es")
plot(inf.analysis, "i2")


#-----------------------------------------------------------------------------#
# functional annotation subset ---- 
# focus on pooled scores for now
#-----------------------------------------------------------------------------#

svm_pooled <- readRDS("/media/work/svm_scores/svm_pooled_filter_sd3.rds")
x1 <- gxe %>%
  filter(SNP %in% svm_pooled$SNP)


# output bin SNPs for expectation based hybrid method..  
twostep_eh_snps(x1, 'chiSqG', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
twostep_eh_snps(x1, 'chiSqGE', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))
twostep_eh_snps(x1, 'chiSqEDGE', output_dir = glue('/media/work/gwis/twostep_expectation_hybrid/{exposure}/svm_subset/'))


# create manhattan, qq, two-step
plot_funcannot_wrap(gxe, exposure = exposure, covariates = covariates, output_dir = output_dir, filename_suffix = "_functional_subset")

# run expectation based hybrid?
simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqG', output_dir = output_dir, filename_suffix = "_functional_subset")
simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqGE', output_dir = output_dir, filename_suffix = "_functional_subset")
simplem_wrap2(x = x1, exposure = exposure, covariates = covariates, simplem_step1_statistic = 'chiSqEDGE', output_dir = output_dir, filename_suffix = "_functional_subset")




#-----------------------------------------------------------------------------#
# GxE additional analysis ---- 
#-----------------------------------------------------------------------------#



# interaction model stratified by study design, sex, and tumor site (NOT 3 way interaction)
snps <- c("6:12577203:T:C", "5:40252294:C:T")
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'study_design', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'sex', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'cancer_site_sum2', path = glue("{path}/output")))

covariates_multi = c(covariates, "smk_ever", "bmi", "alcoholc", "redmeatqc2")
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates_multi, method = "chiSqGxE", strata = 'study_design', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates_multi, method = "chiSqGxE", strata = 'sex', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates_multi, method = "chiSqGxE", strata = 'cancer_site_sum2', path = glue("{path}/output")))






# Additional adjustment covariates
covariates_sets <- list(covariates,
                        c(covariates, 'bmi', 'smk_ever', 'alcoholc', 'redmeatqc2'))

fit_gxe_covars(data = input_data, exposure = exposure, snp = "6:12577203:T:C", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))
fit_gxe_covars(data = input_data, exposure = exposure, snp = "5:40252294:C:T", covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output"))



wtf <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr6_12577203.qs")



# output RERI plots
snps <- c("6:12577203:T:C")
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))

# need to flip exposure
tmp_flip <- input_data %>% 
  mutate(asp_ref = fct_relevel(asp_ref, "Yes"))
snps <- c( "6:32560631:C:T", "5:40252294:C:T")
walk(snps, ~ reri_wrapper(data_epi = tmp_flip, exposure = exposure, snp = .x, covariates = covariates, path = glue("{path}/output")))





# stratified odds ratio table with additional covariates
# covariates <- c(covariates, 'bmi', 'smk_ever', 'alcoholc', 'redmeatqc2')
snps <- c( "6:12577203:T:C", "5:40252294:C:T")
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


covariates_mult <- c(covariates, 'bmi', 'smk_ever', 'alcoholc', 'redmeatqc2')
snps <- c( "6:12577203:T:C", "5:40252294:C:T")
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates_mult,
  path = glue("{path}/output/")))



tmp <- readRDS(glue("/media/work/gwis_test/{exposure}/output/posthoc/stratified_oddsratio_{exposure}_v2.4_chr6_12577203_T_C_age_ref_imp_pc1_pc2_pc3_sex_study_gxe.rds")) %>% 
  setNames(make.unique(names(.))) %>% 
  separate(N0, into = c("N0_case", "N0_control"), sep = "/") %>% 
  separate(N1, into = c("N1_case", "N1_control"), sep = "/") %>% 
  separate(N2, into = c("N2_case", "N2_control"), sep = "/") %>% 
  mutate(across(.cols = starts_with("N"), ~ round(as.numeric(.x), 0) ))


tmp <- readRDS("/media/work/gwis_test/asp_ref/output/posthoc/stratified_oddsratio_asp_ref_v2.4_chr6_12577203_T_C_age_ref_imp_alcoholc_bmi_pc1_pc2_pc3_redmeatqc2_sex_smk_ever_study_gxe.rds") %>% 
  setNames(make.unique(names(.))) %>% 
  separate(N0, into = c("N0_case", "N0_control"), sep = "/") %>% 
  separate(N1, into = c("N1_case", "N1_control"), sep = "/") %>% 
  separate(N2, into = c("N2_case", "N2_control"), sep = "/") %>% 
  mutate(across(.cols = starts_with("N"), ~ round(as.numeric(.x), 0) ))




# SNP information
snp_info <- qread("/media/work/FIGI_RsqEstimate_chrALL.qs") %>% 
  filter(id %in%  c("6:12577203:T:C", "6:32560631:C:T", "5:40252294:C:T"))

ff <- do.call(c, lapply(strsplit(snps, split = ":"), function(x) paste(x[1], as.numeric(x[2])-1, x[2],  sep = ':')))

snp_mart = useMart(biomart = "ENSEMBL_MART_SNP", 
                   host    = "grch37.ensembl.org", 
                   path    = "/biomart/martservice", 
                   dataset = "hsapiens_snp")

rsid = getBM(attributes = c("refsnp_id", "allele", "chr_name", "chrom_end"),
             filters = c("chromosomal_region"),
             values = ff, mart=snp_mart)


snp_info$rsid <- rsid$refsnp_id

snp_info_out <- snp_info %>% 
  separate(SNP, into = c("chr", "bp", "REF", "ALT"), remove= F) %>%
  rename(ALT_AF = GxESet_AltAlleleFreq, 
         MAF = maf, 
         Imputation_Rsq = GxESet_Rsq) %>% 
  dplyr::select(SNP, rsid, REF, ALT, ALT_AF , MAF, Imputation_Rsq )

saveRDS(snp_info_out, glue("{path}/output/posthoc/gwis_snp_info.rds"))



# need allele frequency by study_gxe to confirm finding (and sensitivity analysis)
chr6_125 <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr6_12577203.qs") %>% 
  inner_join(input_data, 'vcfid')

chr6_125 %>% 
  summarise(total = n(), 
            study_aaf = sum(chr6_12577203_T_C_dose) / (total*2))

aaf <- function(x) {
  sum(x) / nrow(x)
}

out <- chr6_125 %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr6_12577203_T_C_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")

ggsave(filename = "~/Dropbox/asp_ref_chr6_12577203_T_C_AAF.png", width = 7.5, height = 4.72)



# would results change if you remove mecc (note wald statistic)
model1 <- glm(glue("outcome ~ {exposure}*chr6_12577203_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_125, family = 'binomial')
summary(model1)

chr6_125b <- dplyr::filter(chr6_125, !grepl("MECC", study_gxe))
model2 <- glm(glue("outcome ~ {exposure}*chr6_12577203_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_125b, family = 'binomial')
summary(model2)


model3 <- glm(glue("outcome ~ chr6_12577203_T_C_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_125, family = 'binomial')
summary(model3)







# need allele frequency by study_gxe to confirm finding (and sensitivity analysis)
chr6_325 <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr6_32560631.qs") %>% 
  inner_join(input_data, 'vcfid')

aaf <- function(x) {
  sum(x) / nrow(x)
}

out <- chr6_325 %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr6_32560631_C_T_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")

ggsave(filename = "~/Dropbox/asp_ref_chr6_32560631_C_T_AAF.png", width = 7.5, height = 4.72)


# would results change if you remove mecc (note wald statistic)
# EXCLUDE !!!!!!!!!!!!!
model1 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_325, family = 'binomial')
summary(model1)

chr6_325b <- dplyr::filter(chr6_325, !grepl("UKB", study_gxe))
model2 <- glm(glue("outcome ~ {exposure}*chr6_32560631_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr6_325b, family = 'binomial')
summary(model2)







# need allele frequency by study_gxe to confirm finding (and sensitivity analysis)
chr5 <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr5_40252294.qs") %>% 
  inner_join(input_data, 'vcfid')

aaf <- function(x) {
  sum(x) / nrow(x)
}

out <- chr5 %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr5_40252294_C_T_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")

ggsave(filename = "~/Dropbox/asp_ref_chr5_40252294_C_T_AAF.png", width = 7.5, height = 4.72)



# would results change if you remove mecc (note wald statistic)
# EXCLUDE!!!!!!!!!!!!!!
model1 <- glm(glue("outcome ~ {exposure}*chr5_40252294_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr5, family = 'binomial')
summary(model1)

chr5b <- dplyr::filter(chr5, !grepl("REACH", study_gxe))
model2 <- glm(glue("outcome ~ {exposure}*chr5_40252294_C_T_dose + {glue_collapse(covariates, sep = '+')}"), data = chr5b, family = 'binomial')
summary(model2)











# ================================================================== #
# ======= rmarkdown reports ---- 
# ================================================================== #
main_effects_report(exposure = exposure, hrc_version = hrc_version, covariates = covariates, path = path)


gwis_report(exposure = exposure, 
            hrc_version = hrc_version, 
            covariates = covariates)

posthoc_report(exposure = exposure, 
               hrc_version = hrc_version,
               covariates = covariates,
               path = path)








# ================================================================== #
# ======= functional follow? ---- 
# ================================================================== #

hrc <- readRDS("~/data/HRC_rsID/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr5.rds")

# Chromosome 5
tmp1 <- fread(glue("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_{exposure}_{hrc_version}_chr5_40252294.ld"))
tmp2 <- filter(hrc, POS %in% tmp1$BP_B) %>% 
  mutate(V1 = paste0("chr", `#CHROM`), 
         V2 = POS - 1, 
         V3 = POS, 
         V4 = paste(`#CHROM`, POS, REF, ALT, sep = "-")) %>% 
  dplyr::select(V1, V2, V3, V4)

write.table(tmp2, file = glue("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_{exposure}_{hrc_version}_chr5_40252294.bed"), quote = F, row.names = F, col.names = F, sep = '\t')



gtex <- fread("~/Dropbox/FIGI/Annotation_Workflow/example_output/figi_asp_ref_v2.4_chr5_40252294_eQTL_overlap_summary.tsv",  header = F)
vep <- fread("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_asp_ref_v2.4_chr5_40252294_vep.txt",  header = T)
out <- inner_join(gtex, vep, by = c("V1"="#Uploaded_variation")) %>% 
  mutate(BP_B = as.integer(unlist(lapply(strsplit(V1, split = '-'), function(x) x[2])) )) %>% 
  inner_join(tmp1[, c('BP_B', 'R2')], by = "BP_B") %>% 
  arrange(desc(R2))


write_tsv(out, file = "~/Dropbox/FIGI/Annotation_Workflow/example_output/figi_asp_ref_v2.4_chr5_40252294_output.tsv", quote = F)



# Chromosome 6
hrc <- readRDS("~/data/HRC_rsID/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr6.rds")
tmp1 <- fread(glue("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_{exposure}_{hrc_version}_chr6_12577203.ld"))
tmp2 <- filter(hrc, POS %in% tmp1$BP_B) %>% 
  mutate(V1 = paste0("chr", `#CHROM`), 
         V2 = POS - 1, 
         V3 = POS, 
         V4 = paste(`#CHROM`, POS, REF, ALT, sep = "-")) %>% 
  dplyr::select(V1, V2, V3, V4)

write.table(tmp2, file = glue("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_{exposure}_{hrc_version}_chr6_12577203.bed"), quote = F, row.names = F, col.names = F, sep = '\t')



gtex <- fread("~/Dropbox/FIGI/Annotation_Workflow/example_output/figi_aspirin_v2.4_chr6_12577203_eQTL_overlap_summary.tsv",  header = F)
vep <- fread("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_aspirin_v2.4_chr6_12577203_vep.txt",  header = T)
out <- inner_join(gtex, vep, by = c("V1"="#Uploaded_variation"))







# Chromosome 6
hrc <- readRDS("~/data/HRC_rsID/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.chr6.rds")
tmp1 <- fread(glue("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_{exposure}_{hrc_version}_chr6_32560631.ld"))
tmp2 <- filter(hrc, POS %in% tmp1$BP_B) %>% 
  mutate(V1 = paste0("chr", `#CHROM`), 
         V2 = POS - 1, 
         V3 = POS, 
         V4 = paste(`#CHROM`, POS, REF, ALT, sep = "-")) %>% 
  dplyr::select(V1, V2, V3, V4)

write.table(tmp2, file = glue("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_{exposure}_{hrc_version}_chr6_32560631.bed"), quote = F, row.names = F, col.names = F, sep = '\t')



gtex <- fread("~/Dropbox/FIGI/Annotation_Workflow/example_output/figi_aspirin_v2.4_chr6_32560631_eQTL_overlap_summary.tsv",  header = F)
vep <- fread("~/Dropbox/FIGI/Annotation_Workflow/example_input/figi_aspirin_v2.4_chr6_32560631_vep.txt",  header = T)
out <- inner_join(gtex, vep, by = c("V1"="#Uploaded_variation")) %>% 
  mutate(BP_B = as.integer(unlist(lapply(strsplit(V1, split = '-'), function(x) x[2])) )) %>% 
  inner_join(tmp1[, c('BP_B', 'R2')], by = "BP_B") %>% 
  arrange(desc(R2))

write_tsv(out, file = "~/Dropbox/FIGI/Annotation_Workflow/example_output/figi_asp_ref_v2.4_chr6_32560631_output.tsv", quote = F)










# ================================================================== #
# updated PCs and how it changes results ------
# ================================================================== #
new_pc <- fread("/media/work/gwis_test/PCA/20210602/figi_gwas_pca_v2.3_EUR_131677.eigenvec") %>% 
  rename(vcfid = `#FID`)

figi <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr6_12577203.qs") %>% 
  inner_join(input_data, 'vcfid')

figi <- qread("/media/work/gwis_test/asp_ref/output/posthoc/dosage_chr5_40252294.qs") %>% 
  inner_join(input_data, 'vcfid')

# just drops one individual, we'll live. 
figi_newpc <- figi %>% 
  dplyr::select(-starts_with("pc")) %>% 
  inner_join(new_pc, 'vcfid')




# original GxE findings (chr6)
modelb <- glm(outcome ~ chr6_12577203_T_C_dose+asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
model1 <- glm(outcome ~ chr6_12577203_T_C_dose*asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
lrtest(modelb, model1)

modelb <- glm(outcome ~ chr6_12577203_T_C_dose+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr6_12577203_T_C_dose*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
lrtest(modelb, model1)


# chr5_40252294_C_T finding (two step expectation based.. not a 1:1 comparison but just to get an idea)

# marginal association
modelb <- glm(outcome ~                          asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr5_40252294_C_T_dose + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
lrtest(modelb, model1)
# from p-value 1.0930e-06 to 1.085e-06 using new PCs


# GxE
modelb <- glm(outcome ~ chr5_40252294_C_T_dose+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr5_40252294_C_T_dose*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
lrtest(modelb, model1)
# from p-value 4.4147e-05 to 4.594e-05 using new PCs 


