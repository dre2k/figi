test <- figi %>% 
  filter(grepl("CCFR_4", study_gxe))
table(test$study_gxe)
table(test$calcium_totqc2)


m1 <- glm(outcome ~ calcium_totqc2 + age_ref_imp + sex + energytot_imp, data = test, family = 'binomial')
exp(coef(m1))





test2 <- figi %>% 
  filter(grepl("CCFR_1", study_gxe))
table(test$study_gxe)

m2 <- glm(outcome ~ calcium_totqc2 + age_ref_imp + sex + energytot_imp, data = test2, family = 'binomial')
exp(coef(m2))


# use original data to make sure

idlist <- test$pooledcompassid

ccfr <- fread("~/data/FIGI_samplefile_epi-200309/epi/101ccfr0.csv") %>% 
  filter(pooledcompassid %in% idlist) %>% 
  mutate(outcome = ifelse(outc == "Control", 0, 1))
table(ccfr$outcome)
table(ccfr$calcium_totqc2)

m3 <- glm(outcome ~ calcium_totqc2 + age_ref + sex + energytot, data = ccfr, family = 'binomial')
exp(coef(m3))






#------------------- check original data vs current ------------------------- #
figi_ccfr4 <- figi %>% 
  filter(grepl("CCFR_4", study_gxe))

table(figi_ccfr4$study_gxe)
table(figi_ccfr4$outcome)
table(figi_ccfr4$calcium_totqc2)
table(figi_ccfr4$sex)


original_ccfr4 <- fread("~/data/FIGI_samplefile_epi-200309/epi/101ccfr0.csv") %>% 
  filter(pooledcompassid %in% idlist) %>% 
  mutate(calcium_totqc2 = calcium_totqc2 - 1, 
         outcome = ifelse(outc == "Control", 0, 1))

table(original_ccfr4$outcome)
table(original_ccfr4$calcium_totqc2)
table(original_ccfr4$sex)

m1 <- glm(outcome ~ calcium_totqc2 + age_ref_imp + sex + energytot_imp, data = figi_ccfr4, family = 'binomial')
exp(coef(m1))

m12 <- glm(outcome ~ calcium_totqc2 + age_ref_imp + sex + energytot, data = figi_ccfr4, family = 'binomial')
exp(coef(m12))

m2 <- glm(outcome ~ calcium_totqc2 + age_ref + sex + energytot, data = original_ccfr4, family = 'binomial')
exp(coef(m2))



# energy might be the problem..
summary(figi_ccfr4$energytot_imp)
summary(figi_ccfr4$energytot)

summary(original_ccfr4$energytot)
summary(figi_ccfr4$energytot)

overlap <- inner_join(figi_ccfr4, original_ccfr4, 'pooledcompassid')

plot(overlap$energytot_imp, overlap$energytot.x)
boxplot(overlap$energytot_imp, overlap$energytot.x)

# inclusion of the ~ 1200 energytot NAs drive ORs down

overlap <- overlap %>% 
  mutate(energy_available = ifelse(is.na(energytot.x), 'no', 'yes'))
table(overlap$energy_available)
table(overlap$calcium_totqc2.x, overlap$energy_available)


ggplot(data=overlap, aes(x=energy_available, fill=as.character(calcium_totqc2.x))) +
  geom_bar(stat="count", position=position_dodge())

# cases/controls by availability of energy..
ggplot(data=overlap, aes(x=energy_available, fill=outc.x)) +
  geom_bar(stat="count", position=position_dodge())
