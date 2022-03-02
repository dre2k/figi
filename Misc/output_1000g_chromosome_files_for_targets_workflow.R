# AK
# mid october, 2021
#
# for functional annotation, one portion was using HRC obtaining ref/alt alleles (for VEP tool), but I think i should stick to 1000G throughout for consistency
# create chromosome specific files that you can use for your workflow

library(data.table)
library(tidyverse)

x <- fread("/media/work/ReferencePanels/1000GP_Phase3_combined.legend")

y <- dplyr::select(x, chr, position, a0, a1)

saveRDS(y %>% dplyr::filter(chr == 1), file = "/media/work/gwis_test/functional_annotation/kgp/1000GP_Phase3_combined.legend_chr1.rds")

for(xx in 3:22) {
  saveRDS(y %>% dplyr::filter(chr == xx), file = glue("/media/work/gwis_test/functional_annotation/kgp/1000GP_Phase3_combined.legend_chr{xx}.rds"))
  
}