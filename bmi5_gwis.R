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
exposure = 'bmi5'
hrc_version = 'v2.3'
annotation_file <- 'gwas_200_ld_annotation_feb2021.txt'
covariates <- sort(c('age_ref_imp', 'sex', 'study_gxe', 'pc1', 'pc2', 'pc3'))
path = glue("/media/work/gwis_test/{exposure}/")


# input data
esubset <- readRDS(glue("/media/work/gwis_test/{exposure}/data/FIGI_{hrc_version}_gxeset_{exposure}_basic_covars_glm.rds")) %>% 
  pull(vcfid)

input_data <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) %>% 
  filter(vcfid%in% esubset) %>% 
  mutate(bmic4f = case_when(bmi < 18.5 ~ "Underweight", 
                            between(bmi, 18.5, 24.99) ~ "Normal",
                            between(bmi, 25, 29.99) ~ "Overweight", 
                            bmi > 29.99 ~ "Obese"),
         bmic3f = fct_relevel(factor(bmic3), "Normal", "Overweight", "Obese"))



#-----------------------------------------------------------------------------#
# main effects ----
#-----------------------------------------------------------------------------#

# ------ meta-analysis ------ #
output_dir = as.character(glue("{path}/output/posthoc/"))
covariates_meta <- sort(covariates[which(!covariates %in% c(paste0(rep('pc', 20), seq(1, 20)), "study_gxe"))])

create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "all", forest_height = 15, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "proximal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "distal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "rectal", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "female", forest_height = 13, categorical = F)
create_forest_plot(data_epi = input_data, exposure = exposure, covariates = covariates_meta, hrc_version = hrc_version, path = glue("{path}/output/posthoc/"), strata = "male", forest_height = 13, categorical = F)


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

snps <- c("15:33122966:C:T", "15:33120215:T:C")

# output GxE models adjusted by different covariate sets
covariates_sets <- list(covariates, 
                        c(covariates, "smk_ever"))
walk(snps, ~ fit_gxe_covars(data_epi = input_data, exposure = exposure, snp = .x, covariates_list = covariates_sets, method = 'chiSqGxE', path = glue("{path}/output")))






# stratified by tumor and sex 

walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'sex', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'cancer_site_sum2', path = glue("{path}/output")))
walk(snps, ~ fit_gxe_stratified(data_epi = input_data, exposure = exposure, snp = .x, covariates = covariates, method = "chiSqGxE", strata = 'study_design', path = glue("{path}/output")))



# stratified odds ratio table with and without smoking
# covariates_multi <- c(covariates, 'smk_ever')
snps <- c( "15:33122966:C:T", "15:33120215:T:C")
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


covariates_multi <- c(covariates, 'smk_ever')
walk(snps, ~ fit_stratified_or(
  data_epi = input_data,
  exposure = exposure,
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates_multi,
  path = glue("{path}/output/")))





# stratified odds ratio table with and without smoking
# covariates_multi <- c(covariates, 'smk_ever')
snps <- c( "15:33122966:C:T", "15:33120215:T:C")
walk(snps, ~ fit_stratified_or_q3(
  data_epi = input_data,
  exposure = 'bmic3f',
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates,
  path = glue("{path}/output/")))


covariates_multi <- c(covariates, 'smk_ever')
walk(snps, ~ fit_stratified_or_q3(
  data_epi = input_data,
  exposure = 'bmic3f',
  snp = .x,
  hrc_version = hrc_version,
  covariates = covariates_multi,
  path = glue("{path}/output/")))







# output RERI plots
snps <- c("15:33122966:C:T", "15:33120215:T:C")
walk(snps, ~ reri_wrapper(data_epi = input_data, exposure = 'bmi3cf', snp = .x, covariates = covariates, path = glue("{path}/output")))




# SNP information
snp_info <- qread("/media/work/FIGI_RsqEstimate_chrALL.qs") %>% 
  filter(id %in%  c("15:33122966:C:T"))

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
tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr15_33122966.qs")) %>% 
  inner_join(input_data, 'vcfid')

aaf <- function(x) {
  sum(x) / nrow(x)
}

out <- tmp %>% 
  group_by(study_gxe) %>% 
  summarise(total = n(), 
            study_aaf = sum(chr15_33122966_C_T_dose) / (total*2)) %>% 
  arrange(study_aaf) %>% 
  mutate(study_gxe = fct_reorder(study_gxe, study_aaf))

ggplot(aes(x = study_gxe, y = study_aaf), data = out) + 
  geom_point() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 270)) + 
  xlab("Study") + 
  ylab("Alternate Allele Frequency")


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


posthoc_report(exposure = exposure)

# test
rmarkdown::render(glue("/home/rak/git/figi/{exposure}_posthoc_test.Rmd"), 
                  output_file = glue("~/Dropbox/FIGI/Results/{exposure}_posthoc_test.html"))


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
# updated PCs and how it changes results.. 
# ================================================================== #

# original GxE findings.. 
modelb <- glm(outcome ~ chr6_12577203_T_C+asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
model1 <- glm(outcome ~ chr6_12577203_T_C*asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
lrtest(modelb, model1)



# let's add updated PCs .. (keep in mind they're calculated using HRC v3.0)
pcs <- fread("~/figi_gxe_pca_update.eigenvec")
figi_newpc <- inner_join(figi, pcs, by = c('vcfid' = 'IID'))


modelb <- glm(outcome ~ chr6_12577203_T_C+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr6_12577203_T_C*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
summary(model1)

lrtest(modelb, model1)

modelb <- glm(outcome ~ chr6_12577203_T_C+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr6_12577203_T_C*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+ study_gxe, data = figi_newpc, family = 'binomial')
summary(model1)

lrtest(modelb, model1)




#  how about the other hit..
# original GxE findings.. 
  modelb <- glm(outcome ~ chr6_32560631_C_T+asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
model1 <- glm(outcome ~ chr6_32560631_C_T*asp_ref + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi, family = 'binomial')
lrtest(modelb, model1)

# let's add updated PCs .. (keep in mind they're calculated using HRC v3.0)
modelb <- glm(outcome ~ chr6_32560631_C_T+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr6_32560631_C_T*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
summary(model1)

lrtest(modelb, model1)




# chr5_40252294_C_T finding (two step expectation based.. not a 1:1 comparison but just to get an idea)
wtf <- filter(gxe , Location == 40252294) # EDGE chi-square = 24.32772

figi$asp_ref2 <- as.numeric(figi$asp_ref)
figi_newpc$asp_ref2 <- as.numeric(figi_newpc$asp_ref)



# marginal association
modelb <- glm(outcome ~                     asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr5_40252294_C_T + asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
lrtest(modelb, model1)
# 24.275 


modelb <- lm(asp_ref2 ~                     age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi_newpc)
model1 <- lm(asp_ref2 ~ chr5_40252294_C_T + age_ref_imp + sex + pc1 + pc2 + pc3 + study_gxe, data = figi_newpc)
lrtest(modelb, model1)
# 0.5652

pchisq(24.8402, lower.tail = F, df = 2)




# GxE
modelb <- glm(outcome ~ chr5_40252294_C_T+asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
model1 <- glm(outcome ~ chr5_40252294_C_T*asp_ref + age_ref_imp + sex + PC1 + PC2 + PC3 + study_gxe, data = figi_newpc, family = 'binomial')
lrtest(modelb, model1)
