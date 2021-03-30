#-----------------------------------------------------------------------------#
# QQ Plots ----
# (using 'ramwas' package')
#-----------------------------------------------------------------------------#
message("create QQ plots")
create_qqplot_ramwas(x = gxe, exposure = exposure, stat = "chiSqG",       df= 1, filename_suffix = "", output_dir = output_dir)
create_qqplot_ramwas(x = gxe, exposure = exposure, stat = "chiSqGxE",     df= 1, filename_suffix = "", output_dir = output_dir)
create_qqplot_ramwas(x = gxe, exposure = exposure, stat = "chiSq2df",     df= 2, filename_suffix = "", output_dir = output_dir)
create_qqplot_ramwas(x = gxe, exposure = exposure, stat = "chiSq3df",     df= 3, filename_suffix = "", output_dir = output_dir)
create_qqplot_ramwas(x = gxe, exposure = exposure, stat = "chiSqGE",      df= 1, filename_suffix = "", output_dir = output_dir)
create_qqplot_ramwas(x = gxe, exposure = exposure, stat = "chiSqControl", df= 1, filename_suffix = "", output_dir = output_dir)
create_qqplot_ramwas(x = gxe, exposure = exposure, stat = "chiSqCase",    df= 1, filename_suffix = "", output_dir = output_dir)



# --------------------------------------------------------------------------- #
# Create manhattan plots
# with and without excluding GWAS, D|G, and E|G (controls) hits (for 2DF/3DF)
#
# edited to include D|G hits (from subset, NOT gwas) in 2df/3df
# --------------------------------------------------------------------------- #


# ------ manhattan plots for each statistic -------- #
# Marginal G Results
create_manhattanplot(x = gxe, exposure, stat = 'chiSqG', output_dir = output_dir, annotation_file = annotation_file)
# GxE results
create_manhattanplot(x = gxe, exposure, stat = 'chiSqGxE', output_dir = output_dir, annotation_file = annotation_file)
# 2DF results
create_manhattanplot(x = gxe, exposure, stat = 'chiSq2df', output_dir = output_dir, annotation_file = annotation_file)
# 3DF results
create_manhattanplot(x = gxe, exposure, stat = 'chiSq3df', output_dir = output_dir, annotation_file = annotation_file)
# GE
create_manhattanplot(x = gxe, exposure, stat = 'chiSqGE', output_dir = output_dir, annotation_file = annotation_file)
# GE Among Controls
create_manhattanplot(x = gxe, exposure, stat = 'chiSqControl', output_dir = output_dir, annotation_file = annotation_file)
# GE Case Only
create_manhattanplot(x = gxe, exposure, stat = 'chiSqCase', output_dir = output_dir, annotation_file = annotation_file)


# ------ exclude gwas -------- #
# these objects should be defined in the parent script:
# - exclude_gwas
# - exclude_g
# - exclude_eg_control
# 
# only going to do exclude_gwas now
gxe_nogwas <- gxe %>%
  dplyr::filter(!SNP2 %in% exclude_gwas)

# # Marginal G Results
# create_manhattanplot(x = gxe_nogwas, exposure, stat = 'chiSqG', output_dir = output_dir, annotation_file = annotation_file, filename_suffix = "_no_gwas")
# # 2DF results
# create_manhattanplot(x = gxe_nogwas, exposure, stat = 'chiSq2df', output_dir = output_dir, annotation_file = annotation_file, filename_suffix = "_no_gwas")
# # 3DF results
# create_manhattanplot(x = gxe_nogwas, exposure, stat = 'chiSq3df', output_dir = output_dir, annotation_file = annotation_file, filename_suffix = "_no_gwas")
# rm(gxe_nogwas)
# gc() 

# this is a quick solution to generate graphs repeatedly and quickly
create_manhattanplot_ramwas(gxe_nogwas, exposure, stat = 'chiSqG', output_dir = output_dir, filename_suffix = "_no_gwas")
create_manhattanplot_ramwas(gxe_nogwas, exposure, stat = 'chiSq2df', output_dir = output_dir, filename_suffix = "_no_gwas")
create_manhattanplot_ramwas(gxe_nogwas, exposure, stat = 'chiSq3df', output_dir = output_dir, filename_suffix = "_no_gwas")
