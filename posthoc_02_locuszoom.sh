#!/bin/bash

#-----------------------------------------------------------------------------#
# locuszoom standalone (CLI)
# general script to run locuszoom on exposures
#
# LD calculation
# - use ASN 1000G_March2012 to calculate LD using axiom_mecc_cfr_ky
# - use AFR 1000G_March2012 to calculate LD using corect_oncoarray CONTROLS 
#   (large file, takes a long time)
# - use EUR to actually use 1000G to calculate LD
# - use ASN 1000G_Nov2014 to use random sample of 1000 FIGI controls 
#   to calculate LD (same one used to clump results for plotting)
#
# for two-step method findings, let's plot chiSqGxE statistic
# (since it should be at least somewhat elevated)
#-----------------------------------------------------------------------------#

# simple modification so you can call this script from R
exposure=${1}
hrc_version=${2}
snp=${3}
stat=${4}

# input p values for plotting
results=/media/work/gwis/locuszoom/FIGI_${hrc_version}_gxeset_${exposure}_basic_covars_gxescan_${stat}_locuszoom.txt

# command call
# using random sample of 1000 FIGI controls
mkdir -p /media/work/gwis/locuszoom/${exposure}
/home/rak/locuszoom/bin/locuszoom --snpset NULL --metal ${results} --pop ASN --build hg19 --source 1000G_Nov2014 --refsnp ${snp} --flank 500kb --prefix /media/work/gwis/locuszoom/${exposure}/locuszoom_${stat}  --no-date title="${exposure} x ${snp} - ${stat}" signifLine=7.30103 signifLineColor="blue" --denote-markers-file /media/work/gwis/locuszoom/locuszoom_annotation_chr${snp}_${stat}.txt


# you just need to make sure that the denote-markers-file has at least ONE marker that's within the regiong being plotted... 
# there are two ways to do this - you can include the marker being plotted... which might actually work out better. 
