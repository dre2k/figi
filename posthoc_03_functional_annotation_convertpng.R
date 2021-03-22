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


#-----------------------------------------------------------------------------#
# assemble .ini file for this SNP, call pygenometracks
# need to output file to the functional_annotation folder
#-----------------------------------------------------------------------------#

filename <-  paste0("/media/work/gwis/functional_annotation/functional_annotation", "_", stat, "_", snp_helper_hit)

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

