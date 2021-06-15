install.packages("interactionR")
library(interactionR)
data (OCdata) ## Case-control data from Rothman and Keller (1972) evaluating the joint effect of alcohol and smoking on oral cancer risk is included in the package (cited in Hosmer and Lemeshow (1992) and Zou (2008))

## fit the interaction model
model.glm <- glm(oc ~ alc*smk, family = binomial(link = "logit"), data = OCdata)

## Then pass the fitted model to the function which generates a list object of class 'interactionR'
value = interactionR(model.glm, exposure_names = c("alc", "smk"), ci.type = "delta", ci.level = 0.95, em = F, recode = F)

interactionR_table(value)



x <- readRDS("~/data/results/input/FIGI_v2.3_gxeset_fruitqc2_basic_covars_glm.rds")

model.glm <- glm(outcome ~ age_ref_imp*sex, family = 'binomial', data = x)
value = interactionR(model.glm, exposure_name = c("age_ref_imp", "sex"), ci.type = "delta", ci.level = 0.05, em = F, recode = T)
interactionR_table(value)





# gwas 140 for comparison with Yi's 

logitP(D=1|G=g,E=e)=β0+β1G+β2E+β3G∗E+β4study_gxe+β5age_ref_imp+β6sex

+β7PC1+β8PC2+β9PC3+β10scenergytot+β11energymiss (if dietary E)



library(tidyverse)
library(data.table)
library(figifs)



12:133701119:G:A
17:10707241:G:A (rs1078643)
5:40280076:G:A (rs12514517)
74280012T>A (rs7121958)
asp_ref <- readRDS("~/data/results/input/FIGI_v2.3_gxeset_asp_ref_basic_covars_glm.rds")[,'vcfid']
exposure_subset <- readRDS(paste0("/home/rak/data/results/input/FIGI_", 'v2.3' , "_gxeset_", 'fiberqc2', "_basic_covars_glm.rds"))[,'vcfid']
x <- readRDS("~/data/results/input/FIGI_v2.3_gxeset_analysis_data_glm.rds") %>% 
  mutate(fiberqc2 = abs(3-fiberqc2)) %>% 
  filter(vcfid %in% exposure_subset)
         # !study_gxe %in% c("HawaiiCCS_AD", "NHS_4", "NHS_5_AD", "REACH_AD", "SMS_AD"))

sort(unique(x$study_gxe))
table(x$outcome, x$study_gxe)



gwas <- readRDS("/media/work/gwis/gwas_hits_140/GWAS_hits_indices_chr11_out.rds")
out <- inner_join(x, gwas, 'vcfid')



studies <- c('ATBC','CCFR_1','CCFR_3', 'CCFR_4', 'CLUEII', 'Colo23', 'CPSII_1', 'CPSII_2', 'CRCGEN', 'DACHS_1', 'DACHS_2', 'DACHS_3', 'DACHS_4', 'DALS_1', 'DALS_2', 'EDRN', 'HawaiiCCS_AD', 'HPFS_1_2', 'HPFS_3_AD', 'HPFS_4', 'HPFS_5_AD', 'Kentucky', 'LCCS', 'MCCS_1', 'MCCS_2', 'MEC_1', 'MEC_2', 'MECC_1', 'MECC_2', 'MECC_3', 'NCCCSI', 'NCCCSII', 'NFCCR_2', 'NHS_1_2', 'NHS_3_AD', 'NHS_4', 'NHS_5_AD', 'PHS', 'PLCO_1_Rematch', 'PLCO_2', 'PLCO_3', 'PLCO_4_AD' , 'REACH_AD', 'SELECT' ,'SMS_AD', 'UKB_1' ,'USC_HRT_CRC', 'VITAL', 'WHI_1', 'WHI_2', 'WHI_3')







wtf <- glm(outcome ~ fiberqc2 + sex + age_ref_imp + study_gxe, data = out, family = 'binomial')
summary(wtf)

wtf <- glm(outcome ~ X11.74280012 + study_gxe + sex + age_ref_imp + pc1 + pc2 + pc3 , data = out, family = 'binomial')
summary(wtf)

wtf <- glm(outcome ~ fiberqc2 * X11.74280012 + study_gxe + sex + age_ref_imp + pc1 + pc2 + pc3 , data = out, family = 'binomial')
summary(wtf)
exp(coef(wtf))

model.glm <- glm(outcome ~ fiberqc2*X11.74280012 + age_ref_imp + study_gxe + pc1 + pc2 + pc3, data = out, family = 'binomial')
value = interactionR(model.glm, exposure_name = c("eo_ref_pm", "X11.74280012"), ci.type = "delta", ci.level = 0.05, em = F, recode = T)
interactionR_table(value)















model <- model.glm

coef(model)

coef = c(2,3,50)
param = 'product'
type = "RERI"
conf.level = 0.95

epi.interaction <- function(model, coef, param = c("product", "dummy"), type = c("RERI", "APAB", "S"), conf.level = 0.95){
  
  N. <- 1 - ((1 - conf.level)/2)
  z <- qnorm(N., mean = 0, sd = 1)
  
  if (class(model)[1] != "glm" & class(model)[2] != "lm" & class(model)[1] != "clogit" & class(model)[1] != "coxph")
    stop("Error: model must be either a glm or coxph object")      
  
  theta1 <- as.numeric(model$coefficients[coef[1]])
  theta2 <- as.numeric(model$coefficients[coef[2]])
  theta3 <- as.numeric(model$coefficients[coef[3]])
  
  if(class(model)[1] == "glm" & class(model)[2] == "lm"){
    theta1.se <- summary(model)$coefficients[coef[1],2]
    theta2.se <- summary(model)$coefficients[coef[2],2]
    theta3.se <- summary(model)$coefficients[coef[3],2]
  }
  
  if(class(model)[1] == "clogit" | class(model)[1] == "coxph"){
    theta1.se <- summary(model)$coefficients[coef[1],3]
    theta2.se <- summary(model)$coefficients[coef[2],3]
    theta3.se <- summary(model)$coefficients[coef[3],3]
  }
  
  if(type == "RERI"){
    
    if(param == "product"){
      cov.mat <- vcov(model)
      
      h1 <- exp(theta1 + theta2 + theta3) - exp(theta1)
      h2 <- exp(theta1 + theta2 + theta3) - exp(theta2)
      h3 <- exp(theta1 + theta2 + theta3)
      
      reri.var <- (h1^2 * theta1.se^2) + 
        (h2^2 * theta2.se^2) + 
        (h3^2 * theta3.se^2) + 
        (2 * h1 * h2 * cov.mat[coef[1],coef[2]]) + 
        (2 * h1 * h3 * cov.mat[coef[1],coef[3]]) + 
        (2 * h2 * h3 * cov.mat[coef[2],coef[3]])
      reri.se <- sqrt(reri.var)
      
      reri.p <- exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1
      reri.l <- reri.p - (z * reri.se)
      reri.u <- reri.p + (z * reri.se)
      rval <- data.frame(est = reri.p, lower = reri.l, upper = reri.u)
    }
    
    if(param == "dummy"){
      cov.mat <- vcov(model)
      
      h1 <- -exp(theta1)
      h2 <- -exp(theta2)
      h3 <-  exp(theta3)
      
      reri.var <- (h1^2 * (cov.mat[coef[1],coef[1]])) + 
        (h2^2 * (cov.mat[coef[2],coef[2]])) + 
        (h3^2 * (cov.mat[coef[3],coef[3]])) + 
        (2 * h1 * h2 * cov.mat[coef[1],coef[2]]) + 
        (2 * h1 * h3 * cov.mat[coef[1],coef[3]]) + 
        (2 * h2 * h3 * cov.mat[coef[2],coef[3]])
      reri.se <- sqrt(reri.var)
      
      reri.p <- exp(theta3) - exp(theta1) - exp(theta2) + 1
      reri.l <- reri.p - (z * reri.se)
      reri.u <- reri.p + (z * reri.se)
      
      rval <- data.frame(est = reri.p, lower = reri.l, upper = reri.u)
    }
  }
  
  if(type == "APAB"){
    
    if(param == "product"){
      cov.mat <- vcov(model)
      
      h1 <- ((exp(theta1 + theta2 + theta3) - exp(theta1)) / (exp(theta1 + theta2 + theta3))) - ((exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1) / (exp(theta1 + theta2 + theta3)))
      h2 <- ((exp(theta1 + theta2 + theta3) - exp(theta2)) / (exp(theta1 + theta2 + theta3))) - ((exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1) / (exp(theta1 + theta2 + theta3)))
      h3 <- 1 -((exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1) / exp(theta1 + theta2 + theta3))
      
      apab.var <- (h1^2 * theta1.se^2) + 
        (h2^2 * theta2.se^2) + 
        (h3^2 * theta3.se^2) + 
        (2 * h1 * h2 * cov.mat[coef[1],coef[2]]) + 
        (2 * h1 * h3 * cov.mat[coef[1],coef[3]]) + 
        (2 * h2 * h3 * cov.mat[coef[2],coef[3]])
      apab.se <- sqrt(apab.var)
      
      apab.p <- (exp(theta1 + theta2 + theta3) - exp(theta1) - exp(theta2) + 1) / exp(theta1 + theta2 + theta3)
      
      apab.l <- apab.p - (z * apab.se)
      apab.u <- apab.p + (z * apab.se)
      rval <- data.frame(est = apab.p, lower = apab.l, upper = apab.u)
    }
    
    if(param == "dummy"){
      cov.mat <- vcov(model)
      
      h1 <- -exp(theta1 - theta3)
      h2 <- -exp(theta2 - theta3)
      h3 <- (exp(theta1) + exp(theta2) - 1) / exp(theta3)
      
      apab.var <- (h1^2 * (cov.mat[coef[1],coef[1]])) + 
        (h2^2 * (cov.mat[coef[2],coef[2]])) + 
        (h3^2 * (cov.mat[coef[3],coef[3]])) + 
        (2 * h1 * h2 * cov.mat[coef[1],coef[2]]) + 
        (2 * h1 * h3 * cov.mat[coef[1],coef[3]]) + 
        (2 * h2 * h3 * cov.mat[coef[2],coef[3]])
      apab.se <- sqrt(apab.var)
      
      # apab.p <- exp(-theta3) - exp(theta1 - theta3) - exp(theta2 - theta3) + 1
      # Equation 4 (Skrondal 2003):
      apab.p <- (exp(theta3) - exp(theta1) - exp(theta2) + 1) / exp(theta3)
      apab.l <- apab.p - (z * apab.se)
      apab.u <- apab.p + (z * apab.se)
      rval <- data.frame(est = apab.p, lower = apab.l, upper = apab.u)
    }
  }
  
  if(type == "S"){
    
    if(param == "product"){
      # Calculate s.p:
      s.p <- (exp(theta1 + theta2 + theta3) - 1) / (exp(theta1) + exp(theta2) - 2)
      cov.mat <- vcov(model)
      
      # If model type is glm or cph and point estimate of S is negative terminate analysis.
      # Advise user to use a linear odds model:
      if(class(model)[1] == "glm" & class(model)[2] == "lm" & s.p < 0){
        message <- paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
        stop(message)
      }
      
      if(class(model)[1] == "clogit" & class(model)[2] == "coxph" & s.p < 0){
        message <- paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
        stop(message)
      }
      
      h1 <- ((exp(theta1 + theta2 + theta3)) / (exp(theta1 + theta2 + theta3) - 1)) - (exp(theta1) / (exp(theta1) + exp(theta2) - 2))
      h2 <- ((exp(theta1 + theta2 + theta3)) / (exp(theta1 + theta2 + theta3) - 1)) - (exp(theta2) / (exp(theta1) + exp(theta2) - 2))
      h3 <- exp(theta1 + theta2 + theta3) / (exp(theta1 + theta2 + theta3) - 1)
      
      lns.var <- h1^2 * theta1.se^2 + 
        h2^2 * theta2.se^2 + 
        h3^2 * theta3.se^2 + 
        (2 * h1 * h2 * cov.mat[coef[2],coef[1]]) + 
        (2 * h1 * h3 * cov.mat[coef[3],coef[1]]) + 
        (2 * h2 * h3 * cov.mat[coef[3],coef[2]])
      
      lns.se <- sqrt(lns.var)
      
      lns.p <- log(s.p)
      lns.l <- lns.p - (z * lns.se)
      lns.u <- lns.p + (z * lns.se)
      
      s.l <- exp(lns.l)
      s.u <- exp(lns.u)
      rval <- data.frame(est = s.p, lower = s.l, upper = s.u)
    }
    
    if(param == "dummy"){
      # Calculate s.p:
      s.p <- (exp(theta3) - 1) / (exp(theta1) + exp(theta2) - 2)
      cov.mat <- vcov(model)
      
      # If model type is glm or cph and point estimate of S is negative terminate analysis.
      # Advise user to use a linear odds model:
      if(class(model)[1] == "glm" & class(model)[2] == "lm" & s.p < 0){
        message <- paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
        stop(message)
      }
      
      if(class(model)[1] == "clogit" & class(model)[2] == "coxph" & s.p < 0){
        message <- paste("Point estimate of synergy index (S) is less than zero (", round(s.p, digits = 2), ").\n  Confidence intervals cannot be calculated using the delta method. Consider re-parameterising as linear odds model.", sep = "")
        stop(message)
      }
      
      # Use delta method (Hosmer and Lemeshow 1992) if model type is glm, clogit or cph:
      h1 <- -exp(theta1) / (exp(theta1) + exp(theta2) - 2)
      h2 <- -exp(theta2) / (exp(theta1) + exp(theta2) - 2)
      h3 <- exp(theta3) / (exp(theta3) - 1)
      
      lns.var <- h1^2 * theta1.se^2 + 
        h2^2 * theta2.se^2 + 
        h3^2 * theta3.se^2 + 
        (2 * h1 * h2 * cov.mat[coef[2],coef[1]]) + 
        (2 * h1 * h3 * cov.mat[coef[3],coef[1]]) + 
        (2 * h2 * h3 * cov.mat[coef[3],coef[2]])
      
      lns.se <- sqrt(lns.var)
      
      lns.p <- log(s.p)
      lns.l <- lns.p - (z * lns.se)
      lns.u <- lns.p + (z * lns.se)
      
      s.l <- exp(lns.l)
      s.u <- exp(lns.u)
      rval <- data.frame(est = s.p, lower = s.l, upper = s.u)
    } 
  }
  
  return(rval)
}