#=============================================================================#
# helper to create functional annotation figures using pygenometracks
# 
# Comments
# - functional annotation tracks stored in /media/work/gwis/functional_annotation
# - tracks of top hit (vertical line) + SNP in LD (bed track) need to be 
#   recreated for each top hit --> This is the main purpose of this script
# - PRIOTIZE SNPs with R2 > 0.4 (change if collaborators suggest)
# - generate pygenometracks .ini files, then run function call to actually 
#   create the plots
#
# unlike the initial version, i don't loop over SNPs, so generate plots 'manually'
# for each top hit. I think that's cleaner and more purposeful
#=============================================================================#

library(tidyverse)
library(data.table)
# rm(list = ls())

# argument inputs
# exposure = 'alcoholc_moderate'
# snp =  '17:27861663'
# stat = 'chiSq3df'


args <- commandArgs(trailingOnly=T)
exposure <- args[1] # ex: asp_ref
snp <- args[2] # ex: 14:73529409
stat <- args[3] # ex: chiSqGxE


# output directory
# (create this in bash before calling this script..)
dir.create(paste0("/media/work/gwis/functional_annotation_out/", exposure), showWarnings = FALSE)
wdir <-    paste0("/media/work/gwis/functional_annotation_out/", exposure, '/')


#-----------------------------------------------------------------------------#
# read in ld SNPs from locuszoom plot
# (since it generates it for us, might as well use)
#-----------------------------------------------------------------------------#
# !@#&^%&@^%@$
snp_chr <- strsplit(snp, ':')[[1]][1]
snp_bp <- as.numeric(strsplit(snp, ':')[[1]][2])
snp_helper <- paste0('chr', snp_chr, "_", snp_bp-500000)
snp_helper_hit <- paste0('chr', snp_chr, "_", snp_bp)
file_helper <- list.files(paste0('/media/work/gwis/locuszoom/', exposure), recursive = T, full.names = T)
load(file_helper[grepl(snp_helper, file_helper) & grepl(stat, file_helper) & grepl("Rdata",file_helper)])


#-----------------------------------------------------------------------------#
# create bed file + matching pygenometracks .ini call
# create vline file + matching pygenometracks .ini call
#-----------------------------------------------------------------------------#

filename_out <- paste0(wdir, "functional_annotation", "_", stat, "_", snp_helper_hit)

# generate color coded bed track 
# (edited so cutoff is 0.6 rather than 0.4)
# 0.4-0.6 = green
# 0.6-0.8 = red
# > 0.8 = purple
# bed_out <- metal %>% 
#   filter(rsquare >= 0.4) %>% 
#   mutate(v1 = paste0('chr', chr), 
#          v2 = pos_int - 1, 
#          v3 = pos_int,
#          v4 = '.',
#          v5 = round(rsquare*1000), 
#          v6 = '.', 
#          v7 = pos_int - 1, 
#          v8 = pos_int - 1, 
#          v9 = ifelse(rsquare >= 0.4 & rsquare < 0.6, '92,184,92',
#                 ifelse(rsquare >= 0.6 & rsquare < 0.8, '212,63,58', 
#                   ifelse(rsquare >= 0.8, '150,50,184', NA)))) %>% 
#   filter(!duplicated(.)) %>% 
#   dplyr::select(v1, v2,v3,v4,v5,v6,v7,v8, v9)
bed_out <- metal %>% 
  filter(rsquare >= 0.5) %>% 
  mutate(v1 = paste0('chr', chr), 
         v2 = pos_int - 1, 
         v3 = pos_int,
         v4 = paste0(chr, "_", pos_int),
         v5 = round(rsquare*1000), 
         v6 = '.', 
         v7 = pos_int - 1, 
         v8 = pos_int - 1, 
         v9 = ifelse(rsquare >= 0.5 & rsquare < 0.8, '212,63,58', 
                  ifelse(rsquare >= 0.8, '150,50,184', NA))) %>% 
  filter(!duplicated(.)) %>% 
  dplyr::select(v1, v2,v3,v4,v5,v6,v7,v8, v9)



write.table(bed_out, 
            file = paste0(filename_out, '.bed'), 
            quote = F, row.names = F, col.names = F, sep = '\t')

# vlines table
bed_chr_vlines <- bed_out %>% 
  filter(v3 == snp_bp)

write.table(bed_chr_vlines, 
            file = paste0(filename_out, '_vlines.bed'), 
            quote = F, row.names = F, col.names = F, sep = '\t')


# # generate .ini for pygenometracks
# cat(paste("[snps]",
#           paste0("file=", filename_out, '.bed'),
#           "title = r^2 > 0.4",
#           "height = 2",
#           "color = bed_rgb",
#           "border_color=none",
#           "labels=false",
#           "display=collapsed",
#           "file_type=bed", sep = '\n'),
#     file = paste0(filename_out, ".ini"), append = F)

# generate .ini for pygenometracks
#cat(paste("[snps]",
#          paste0("file=", filename_out, '.bed'),
#          "title = r^2 > 0.6",
#          "height = 2",
#          "color = bed_rgb",
#          "border_color=none",
#          "labels=true",
#          "display=stacked",
#          "fontsize=10",
#          "file_type=bed", sep = '\n'),
#    file = paste0(filename_out, ".ini"), append = F)

cat(paste("[snps]",
          paste0("file=", filename_out, '.bed'),
          "title = r^2 > 0.5",
          "height = 1",
          "color = bed_rgb",
          "border_color=none",
          "labels=false",
          "display=collapsed",
          "fontsize=11",
          "file_type=bed", sep = '\n'),
    file = paste0(filename_out, ".ini"), append = F)


cat(paste("\n",
          "[vlines]", 
          paste0("file=", filename_out, '_vlines.bed'),
          "type=vlines", 
          "labels=true", sep = '\n'), # this is only available in updated version of pygenometracks.. 
    file = paste0(filename_out, ".ini"), append = T)



#-----------------------------------------------------------------------------#
# assemble .ini file for this SNP, call pygenometracks
# need to output file to the functional_annotation folder
#-----------------------------------------------------------------------------#

filename <-  paste0("/media/work/gwis/functional_annotation/functional_annotation", "_", stat, "_", snp_helper_hit)

system(paste0("cat /media/work/gwis/functional_annotation/tmp1.ini ", 
              paste0(filename_out, ".ini"), " ",
              "/media/work/gwis/functional_annotation/tmp2.ini ", 
              " > ", paste0(filename, '.ini')))

# create shell script
cat(paste("#!/bin/bash",
          "cd /media/work/gwis/functional_annotation/",
          "\n", 
          paste0("/home/rak/anaconda3/bin/pyGenomeTracks  --tracks ", 
                 paste0(filename, '.ini'), " ", 
                 " --fontSize 8 --dpi 60 ", 
                 " --region chr", snp_chr, ":", snp_bp-50000, "-", snp_bp+50000, 
                 " --outFileName ", 
                 paste0(filename_out, ".pdf")),
          "\n", 
          paste0("gs -o ", 
                 paste0(filename_out, ".png "), " -sDEVICE=png16m -dTextAlphaBits=4 -r300 -dLastPage=1 ", 
                 paste0(filename_out, ".pdf ")), sep = "\n"),
    file = paste0(filename_out, ".sh"), append = F)
    
# call shell script
system(paste0("bash ", filename_out, ".sh"))




#------------------------#
# small final step.. 
# create shell script for converting locuszoom plot to png (purely for convenience)
# ------ > SAVING PNG FILE IN THE FUNCTIONAL ANNOTATION FOLDER!!!!
#------------------------#
file_helper <- list.files(paste0('/media/work/gwis/locuszoom/', exposure), recursive = T, full.names = T)
file_helper_pdf <- (file_helper[grepl(snp_helper, file_helper) & grepl(stat, file_helper) & grepl("pdf",file_helper)])

filename_locuszoom <- paste0(wdir, "locuszoom", "_", stat, "_", snp_helper_hit)

cat(paste("#!/bin/bash",
          "cd /media/work/gwis/functional_annotation/",
          paste0("gs -o ", 
                 paste0(filename_locuszoom, ".png "), " -sDEVICE=png16m -dTextAlphaBits=4 -r300 -dLastPage=1 ", 
                 file_helper_pdf), sep = "\n"),
    file = paste0(filename_out, "_locuszoom.sh"), append = F)

# call shell script
system(paste0("bash ", filename_out, "_locuszoom.sh"))
