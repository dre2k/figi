#!/bin/bash

# females only
#Rscript posthoc_run_analysis.R hrt_ref_pm2 age_ref_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R eo_ref_pm_gxe age_ref_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R ep_ref_pm_gxe age_ref_imp study_gxe pc1 pc2 pc3

# non-dietary (no energytot_imp adjustment)
#Rscript posthoc_run_analysis.R asp_ref       age_ref_imp sex study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R diab          age_ref_imp sex study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R hrt_ref_pm2 age_ref_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R hrt_ref_pm2 age_ref_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R hrt_ref_pm2 age_ref_imp study_gxe pc1 pc2 pc3

# dietary (adjusted by energytot_imp)
#Rscript posthoc_run_analysis.R folate_totqc2 age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R folate_dietqc2 age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R alcoholc_heavy_vs_moderate age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R alcoholc_moderate age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R calcium_totqc2    age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R calcium_dietqc2   age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R fruitqc2          age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R fiberqc2          age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R vegetableqc2      age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R redmeatqc2        age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3
#Rscript posthoc_run_analysis.R procmeatqc        age_ref_imp sex energytot_imp study_gxe pc1 pc2 pc3

#Rscript posthoc_run_analysis.R asp_ref "age_ref_imp+study_gxe+pc1+pc2+pc3+study_gxe:pc1+study_gxe:pc2+study_gxe:pc3+study_gxe:age_ref_imp+study_gxe:p1 + study_gxe:p2+ study_gxe:hrt_ref_pm2"  age_ref_imp study_gxe pc1 pc2 pc3
Rscript posthoc_run_analysis_yi.R asp_ref "asp_ref + p1 + p2 + asp_ref*p1 + asp_ref*p2 + age_ref_imp + study_gxe + pc1 + pc2 + pc3"  age_ref_imp sex study_gxe pc1 pc2 pc3
Rscript posthoc_run_analysis_yi.R hrt_ref_pm2 "hrt_ref_pm2 + p1 + p2 + hrt_ref_pm2*p1 + hrt_ref_pm2*p2 + age_ref_imp + study_gxe + pc1 + pc2 + pc3"  age_ref_imp study_gxe pc1 pc2 pc3
Rscript posthoc_run_analysis_yi.R ep_ref_pm_gxe "ep_ref_pm_gxe + p1 + p2 + ep_ref_pm_gxe:p1 + ep_ref_pm_gxe:p2 + age_ref_imp + study_gxe + pc1 + pc2 + pc3"  age_ref_imp study_gxe pc1 pc2 pc3
Rscript posthoc_run_analysis_yi.R eo_ref_pm_gxe "eo_ref_pm_gxe + p1 + p2 + eo_ref_pm_gxe:p1 + eo_ref_pm_gxe:p2 + age_ref_imp + study_gxe + pc1 + pc2 + pc3"  age_ref_imp study_gxe pc1 pc2 pc3
