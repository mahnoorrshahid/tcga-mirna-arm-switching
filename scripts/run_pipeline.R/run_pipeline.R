
setwd("C:/Users/mahno/OneDrive/Documents/tcga-mirna-arm-switching")

source("scripts/softcode_final.R")

run_arm_bias_pipeline("TCGA_KIRC")
run_all_cancers()