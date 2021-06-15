# check GxSex results after filtering for AAF (CALCULATED EXCLUDING UKB)
# 6/11/2021
# 

gxe <- qread("/media/work/gwis_test/sex/results/FIGI_v2.3_gxeset_sex_basic_covars_gxescan_results.qs")

gxeout <- filter(gxe, SNP %in% aaf_out)

gxe_sig <- filter(gxeout, chiSqGxE_p < 5e-8)



# take GxE clumped 
gxe_clumped <- fread("/media/work/gwis_test/sex/data/FIGI_v2.3_gxeset_sex_chiSqGxE_ldclump.clumped") %>% 
  filter(SNP %in% aaf_out)




figi <- readRDS("/media/work/gwis_test/data/FIGI_v2.3_gxeset_analysis_data_glm.rds") %>% 
  filter(!is.na(sex), 
         # EUR_subset == 1, 
         gxe == 1)

drops <- data.frame(table(figi$studyname, figi$sex)) %>% 
  filter(Freq == 0)

figi_out <- figi %>% 
  filter(!study_gxe %in% unique(drops$Var1))


table(figi$study_gxe, figi$sex)


drops <- data.frame(table(figi$study_gxe, figi$sex)) %>% 
  filter(Freq == 0)

figi_out <- figi %>% 
  filter(!study_gxe %in% unique(drops$Var1))

Var1 Var2 Freq
1         ATBC    0    0
2     HPFS_1_2    0    0
3    HPFS_3_AD    0    0
4       HPFS_4    0    0
5    HPFS_5_AD    0    0
6          PHS    0    0
7       SELECT    0    0
8      NHS_1_2    1    0
9     NHS_3_AD    1    0
10       NHS_4    1    0
11    NHS_5_AD    1    0
12 USC_HRT_CRC    1    0
13       WHI_1    1    0
14       WHI_2    1    0
15       WHI_3    1    0







