hrt <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_hrt_ref_pm2_basic_covars_glm.rds")

x <- filter(input_data, grepl("CCFR|NHS", study_gxe))

table(x$study_gxe, x$hrt_ref_pm2)
table(x$study_gxe, x$eo_ref_pm)
table(x$study_gxe, x$eo_ref_pm_gxe)
table(x$study_gxe, x$ep_ref_pm)
table(x$study_gxe, x$eo_ref_pm_gxe)



xx <- data.frame(table(input_data$eo_ref_pm, input_data$hrt_ref_pm2))
names(xx) = c("eo_ref_pm", "hrt_ref_pm2")
xx

xx <- data.frame(table(input_data$eo_ref_pm_gxe, input_data$hrt_ref_pm2))
names(xx) = c("eo_ref_pm_gxe", "hrt_ref_pm2")
xx


xx <- data.frame(table(input_data$ep_ref_pm, input_data$hrt_ref_pm2))
names(xx) = c("ep_ref_pm", "hrt_ref_pm2")
xx

xx <- data.frame(table(input_data$eo_ref_pm_gxe, input_data$hrt_ref_pm2))
names(xx) = c("eo_ref_pm_gxe", "hrt_ref_pm2")
xx






x <- filter(input_data, grepl("CCFR", study_gxe))

xx <- data.frame(table(x$eo_ref_pm, x$hrt_ref_pm2))
names(xx) = c("eo_ref_pm", "hrt_ref_pm2")
xx

xx <- data.frame(table(x$eo_ref_pm_gxe, x$hrt_ref_pm2))
names(xx) = c("eo_ref_pm_gxe", "hrt_ref_pm2")
xx



xx <- data.frame(table(x$ep_ref_pm, x$hrt_ref_pm2))
names(xx) = c("ep_ref_pm", "hrt_ref_pm2")
xx

xx <- data.frame(table(x$ep_ref_pm_gxe, x$hrt_ref_pm2))
names(xx) = c("ep_ref_pm_gxe", "hrt_ref_pm2")
xx



x <- filter(input_data, grepl("NHS", study_gxe))

xx <- data.frame(table(x$eo_ref_pm, x$hrt_ref_pm2))
names(xx) = c("eo_ref_pm", "hrt_ref_pm2")
xx

xx <- data.frame(table(x$eo_ref_pm_gxe, x$hrt_ref_pm2))
names(xx) = c("eo_ref_pm_gxe", "hrt_ref_pm2")
xx



xx <- data.frame(table(x$ep_ref_pm, x$hrt_ref_pm2))
names(xx) = c("ep_ref_pm", "hrt_ref_pm2")
xx

xx <- data.frame(table(x$ep_ref_pm_gxe, x$hrt_ref_pm2))
names(xx) = c("ep_ref_pm_gxe", "hrt_ref_pm2")
xx



# for NHS, if you noticed, there's no "No's" in the hrt_ref_pm2 variable. they don't have that information. 

eo <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_eo_ref_pm_gxe_basic_covars_glm.rds")
table(eo$study_gxe, eo$eo_ref_pm_gxe)

ep <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_ep_ref_pm_gxe_basic_covars_glm.rds")
table(ep$study_gxe, ep$ep_ref_pm_gxe)




# check rs79439591 for Yu
# 1:53785007:C:T
# for ep_ref_pm_gxe

library(lmtest)
library(expss)

ep <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_ep_ref_pm_gxe_basic_covars_glm.rds")
snps <- readRDS("/media/work/gwis/posthoc/gwis_sig_results_output_ep_ref_pm_gxe.rds")
gxe <- inner_join(ep, snps, 'vcfid')
names(gxe)
snp <- "X1.53785007.C.T_dose"


out <-     glm(glue("outcome ~ ep_ref_pm_gxe*{snp} + age_ref_imp + PC1 + PC2 + PC3 + study_gxe"), data = gxe, family = 'binomial')
out_ref <- glm(glue("outcome ~ ep_ref_pm_gxe+{snp} + age_ref_imp + PC1 + PC2 + PC3 + study_gxe"), data = gxe, family = 'binomial')
lrtest(out, out_ref)


# "pure"?




# first, let's see if we get same estimates as yi..
tmp <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_hrt_ref_pm2_basic_covars_glm.rds") %>% 
  pull(vcfid)

hrt <- readRDS("/media/work/gwis/results/input/FIGI_v2.4_gxeset_analysis_data_glm.rds") %>% 
  filter(!is.na(hrt_ref_pm2))

hrt <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_analysis_data_glm.rds") %>% 
  filter(vcfid %in% tmp)

table(hrt$study, hrt$hrt_ref_pm2) # ccfr small discrepancy wtf





tmp <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_eo_ref_pm_basic_covars_glm.rds") %>% 
  pull(vcfid)

out <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_analysis_data_glm.rds") %>% 
  filter(vcfid %in% tmp)

table(out$study, out$eo_ref_pm) # ccfr small discrepancy wtf




tmp <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_ep_ref_pm_basic_covars_glm.rds") %>% 
  pull(vcfid)

out <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_analysis_data_glm.rds") %>% 
  filter(vcfid %in% tmp)

table(out$study, out$ep_ref_pm) # ccfr small discrepancy wtf



tmp <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_ep_ref_pm_gxe_basic_covars_glm.rds") %>% 
  pull(vcfid)

out <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_analysis_data_glm.rds") %>% 
  filter(vcfid %in% tmp)

table(out$ep_ref_pm_gxe)
table(out$study_gxe, out$ep_ref_pm_gxe) # ccfr small discrepancy wtf




tmp <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_eo_ref_pm_gxe_basic_covars_glm.rds") %>% 
  pull(vcfid)
out <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_analysis_data_glm.rds") %>% 
  filter(vcfid %in% tmp)
table(out$eo_ref_pm_gxe)







figi_gxe <- readRDS("/media/work/gwis/results/input/FIGI_v2.3_gxeset_analysis_data_glm.rds")


# recoding variables.. 


# ep_ref_pm_gxe: 
# control = no EO, no EP, no HRT
# exposed = yes EP, regardless of EO or HRT status

# for instance:
table(figi_gxe$ep_ref_pm)
table(figi_gxe$ep_ref_pm, figi_gxe$hrt_ref_pm2)
ep_yes <- filter(figi_gxe, ep_ref_pm == "Yes")
table(ep_yes$eo_ref_pm, ep_yes$hrt_ref_pm2)\
cro(ep_yes$eo_ref_pm, ep_yes$hrt_ref_pm2)
# No  Yes
# No     0 1955
# Yes    0  493

 # |                  |              | ep_yes$hrt_ref_pm2 |      |
 # |                  |              |                 No |  Yes |
 # | ---------------- | ------------ | ------------------ | ---- |
 # | ep_yes$eo_ref_pm |           No |                    | 1955 |
 # |                  |          Yes |                    |  493 |
 # |                  | #Total cases |                    | 2448 |

ep_no <- filter(figi_gxe, ep_ref_pm == "No")
cro(ep_no$eo_ref_pm, ep_no$hrt_ref_pm2)

 # |                 |              | ep_no$hrt_ref_pm2 |      |
 # |                 |              |                No |  Yes |
 # | --------------- | ------------ | ----------------- | ---- |
 # | ep_no$eo_ref_pm |           No |              6008 |  142 |
 # |                 |          Yes |                   | 2938 |
 # |                 | #Total cases |              6008 | 3080 |
 # 



# 'pure' ep_all:
# control = no EO, no EP, no HRT
# exposed = yes EP, no EO, no HRT
# N should be 1955 exposed, 6008 controls
gxe <- figi_gxe %>% 
  mutate(ep_ref_pm_gxe = ifelse(ep_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                         ifelse(hrt_ref_pm2 == "No", "No", NA)),
         pure_ep_all = ifelse(ep_ref_pm == "Yes" & eo_ref_pm == "No", "Yes",
                       ifelse(ep_ref_pm == "No" & hrt_ref_pm2 == "No" & eo_ref_pm == "No", "No", NA)))

table(gxe$ep_ref_pm_gxe)
table(gxe$pure_ep_all)

table(gxe$study_gxe, gxe$ep_ref_pm_gxe)
table(gxe$study_gxe, gxe$hrt_ref_pm2)
table(gxe$ep_ref_pm_gxe, gxe$pure_ep_all, useNA = 'ifany')








exposure = 'ep_ref_pm_gxe'
min_cell_size = 0
# make sure
gxe_ep_ref_pm_gxe <- figi_gxe %>% 
  filter(!study_gxe %in%  c("HawaiiCCS_AD", "REACH_AD", "SMS_AD"), 
         !is.na(ep_ref_pm)) %>% 
  mutate(ep_ref_pm_gxe = ifelse(ep_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                ifelse(hrt_ref_pm2 == "No", "No", NA)),
         pure_ep_all = ifelse(ep_ref_pm == "Yes" & eo_ref_pm == "No", "Yes",
                              ifelse(ep_ref_pm == "No" & hrt_ref_pm2 == "No" & eo_ref_pm == "No", "No", NA))) 
  # filter(!study_gxe %in% c("CLUEII", "Colo23", "CRCGEN", "DACHS_1", "DACHS_2", "DACHS_3", "DALS_1", "DALS_2", "EPIC", "ESTHER_VERDI", "LCCS", "MCCS_1", "MCCS_2", "MECC_3", "NCCCSII", "NFCCR_2", "NHS_4", "NHS_5", "PLCO_1_Rematch", "PLCO_2", "PLCO_3", "PLCO_4_AD", "SMC_COSM", "UKB_1"))

table(gxe_ep_ref_pm_gxe$study_gxe, gxe_ep_ref_pm_gxe$ep_ref_pm_gxe) 






gxe_ep_ref_pm_gxe <- figi_gxe %>% 
  mutate(ep_ref_pm_gxe = ifelse(ep_ref_pm == "Yes" & hrt_ref_pm2 == "Yes", "Yes", 
                                ifelse(hrt_ref_pm2 == "No", "No", NA)),
         pure_ep_allNo = ifelse(ep_ref_pm == "Yes" & eo_ref_pm == "No", "Yes",
                              ifelse(ep_ref_pm == "No" & hrt_ref_pm2 == "No" & eo_ref_pm == "No", "No", NA))) 

table(gxe_ep_ref_pm_gxe$pure_ep_allNo)
cro(gxe_ep_ref_pm_gxe$pure_ep_allNo)

table(gxe_ep_ref_pm_gxe$study_gxe, gxe_ep_ref_pm_gxe$pure_ep_allNo)





drops <- data.frame(table(gxe_ep_ref_pm_gxe$outcome, gxe_ep_ref_pm_gxe[, exposure], 
                          gxe_ep_ref_pm_gxe$study_gxe)) %>% filter(Freq <= min_cell_size)
tmp <- filter(gxe_ep_ref_pm_gxe, !study_gxe %in% unique(drops$Var3)) %>% 
  dplyr::mutate(study_gxe = fct_drop(study_gxe)) %>% 
  dplyr::select(age_ref_imp, sex, study_gxe, ep_ref_pm) %>% filter(complete.cases(.))

table(tmp$study_gxe, tmp$outcome)
table(tmp$study_gxe, tmp$ep_ref_pm_gxe)


# in summary - your numbers are fine, you need to remove studies with empty cells. 
# just make sure you redefine the ep and eo variable to be even more stringent on who you include 


out <- glm(glue("outcome ~ hrt_ref_pm2 + age_ref_imp + PC1 + PC2 + PC3 + study_gxe"), data = hrt2, family = 'binomial')
summary(out)
