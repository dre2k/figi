library(vcfR)

# x <- read.vcfR("/media/work/gwis_test/misc/output_topmed/rs189056537.vcf")
# out <- vcfR2tidy(x)
# out <- extract_gt_tidy(x , gt_column_prepend = "wtf")


# there were multiallleic markers are manually removed using vim 
# filelist <- list.files("/media/work/gwis_test/misc/output_topmed/", pattern = ".vcf", full.names = T)
filelist_rsid <- gsub(".vcf", "", list.files("/media/work/gwis_test/misc/output_topmed/", pattern = ".vcf"))

x <- lapply(filelist_rsid, function(x) extract_gt_tidy(read.vcfR(glue("/media/work/gwis_test/misc/output_topmed/{x}.vcf")), gt_column_prepend = glue("{x}_")) %>% dplyr::select(Indiv, glue("{x}_DS")))

lapply(x, dim)

y <- Reduce(inner_join, x)

write.csv(y, "~/Dropbox/tdrkh_topmed_dosages.csv", row.names = F)
