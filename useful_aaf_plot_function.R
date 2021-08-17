# function to output MAF plot for convenience

# dat <- input_data
# exposure
# snp = "2:135594399:A:C"

output_aaf_plot <- function(dat, exposure, snp) {
   aaf <- function(x) {
      sum(x) / nrow(x)
   }
   
   # SNP information
   snp_info <- unlist(strsplit(snp, split = ":"))
   
   snpname_clean <- function(x) {
      tmp <- gsub("\\:", "\\_", x)
      # tmp <- gsub("X", "chr", tmp)
      tmp <- glue("chr{tmp}_dose")
      return(tmp)
   }
   
   snpfix <- snpname_clean(snp)
   # snpfix_short <- paste0("chr", gsub("\\:", "\\_", snp))

   tmp <- qread(glue("/media/work/gwis_test/{exposure}/output/posthoc/dosage_chr{snp_info[1]}_{snp_info[2]}.qs")) %>% 
      inner_join(dat, 'vcfid')
   
   maf_compare <- tmp %>% 
      group_by(study_gxe) %>% 
      summarise(total = n(), 
                study_aaf = sum(!! sym(snpfix)) / (total*2)) %>% 
      arrange(study_aaf) %>% 
      mutate(study_gxe = fct_reorder(study_gxe, study_aaf))
   
   ggplot(aes(x = study_gxe, y = study_aaf), data = maf_compare) + 
      geom_point() + 
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 270)) + 
      xlab("Study") + 
      ylab("Alternate Allele Frequency")
   
   ggsave(filename = glue("{path}/output/posthoc/allele_freq_{snpfix}.png"), width = 6, height = 4)
}
