#=============================================================================#
# Nick PRS
# 03/28/2020
#
# Nick calculating PRS for analysis of biomarker data (not FIGI)
# 
# PRS based on CRC don't overlap with data he has available for biomarker data
# so let's find another SNP in close LD
# rs77969132	chr 12 31594813	C	T	-7.9265
#=============================================================================#

figi <- readRDS("~/data/FIGI_BDose_IndexFiles/FIGI_chr12.rds")

head(figi[[11]])

# does this SNP exist in the HRC imputed dataset ...
# yes, the problem is it doesn't exist in the biomarker genotype data
snps <- figi[[11]]
tmp <- filter(snps, Location == 31594813 )

# let's get flanking region 500BP from binarydosage, calculate R2, get a list of SNPs with R2 > 0.9 or so,
# then go to EUR meta-analysis results and output estimates + pvalues etc for him

snp <- "12:31594813:C:T"

out <- snps %>% 
    mutate(id = paste(Chromosome, Location, Reference, Alternate, sep = ":"))
  
tmp <- out %>% 
    filter(Location >= 31594813-100000 & Location <= 31594813+100000)
  
tmp_index <- which(out$id %in% tmp$id)
  
saveRDS(tmp_index, file = paste0("~/nick_chr12_31594813.rds"), version = 2)




#--------------------------------------------------#
# After calculating r2 with SNPs +- 100k bp
#--------------------------------------------------#

snps_r2 <- readRDS("~/nick_12_31594813_r2.rds")
hist(snps_r2)
snps_r2 <- snps_r2[which(snps_r2 > 0.5)]


# none are r2 > 0.5... 
# take a look at jeroen's ld based annotation file to make sure this is the correct SNP

jeroen <- fread("~/Dropbox/FIGI/Documentation/Annotations_GWASHits_20190913/crc_gwas_indep_signals_050719_temp.tsv")[,-c(16,17)] %>% mutate(SNP = paste(CHROM, POSITION, sep = ":")) %>% 
  filter(CHROM == 12)
