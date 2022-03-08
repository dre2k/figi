rm(list = ls())
load("~/data/FIGI_EpiData_rdata/FIGI_Genotype_Epi_190424.RData")




gxe_set <- readRDS(glue("/media/work/gwis_test/data/FIGI_{hrc_version}_gxeset_analysis_data_glm.rds")) 


# GxE set
gxe_set <- figi %>%
  filter(drop == 0 & gxe == 1)

# PCA Results
pc30k <- fread("~/data/PCA/190506/FIGI_GxESet_KGP_pc20_190430.eigenvec", skip = 1,
               col.names = c("FID", "IID", paste0(rep("PC", 20), seq(1,20)))) %>% 
  mutate(vcfid = IID)

df_tmp <- full_join(gxe_set, pc30k, by="vcfid") %>% 
  dplyr::select(vcfid, outc, race_self, study_gxe, study_site, study, platform, paste0(rep("PC", 20), seq(1,20)))


# add 1000G sample info for plotting
kgp_samples <- fread("/home/rak/data/PCA/integrated_call_samples_v3.20130502.ALL.panel.fix", stringsAsFactors = F) %>%
  dplyr::rename(vcfid = sample)

# final data.frame
df <- full_join(df_tmp, kgp_samples, by = "vcfid") %>% 
  mutate(group = factor(replace(super_pop, is.na(super_pop), 'FIGI'),
                        levels=c("FIGI","AFR", "AMR", "EAS", "EUR", "SAS")), 
         eur_subset = ifelse(vcfid %in% c(gxe_set$vcfid, kgp_samples$vcfid), 1, 0)) %>% 
  filter(eur_subset == 1)






ggplot(data = df, aes(x = PC2, y = PC1, color = group)) +
  geom_point(alpha = 0.5) +
  labs(x = "PC2",
       y = "PC1",
       title = paste0("PC1 vs PC2")) +
  scale_colour_manual(values=c("gold", "red", "black", "purple", "green", "royalblue")) +
  theme_classic() +
  theme(legend.key.size = unit(0.15, 'inches'))


