setwd("...../Survival analysis without CR/")
rm(list = ls())

library(PMCMRplus)
library(PMCMR)

library(ggplot2)
options(scipen=200)

set.seed(100)
# ------------------------------------------- model performance ---------------------
load("result.ci.cox.R")
res.cox <- data.frame("cindex" = result.ci$t[, 2], "ibs" = result.ci$t[, 4])
load("result.ci.aft.R")
res.aft <- data.frame("cindex" = result.ci$t[, 2], "ibs" = result.ci$t[, 4])
load("result.ci.rsf.R")
res.rsf <- data.frame("cindex" = result.ci.rsf$t[, 2], "ibs" = result.ci.rsf$t[, 4])
res.mtlr <- read.csv("../data/results.ci.mtlr.csv")
res.ssmtlr <- read.csv("../data/results.ci.ssmtlr.csv")


auc_results <- data.frame("cox" = res.cox$cindex, 
                          "aft" = res.aft$cindex,
                          "rsf" = res.rsf$cindex,
                          "mtlr" = res.mtlr$cindex)
                     
## Export fle
#library(writexl)
#write_xlsx(auc_results, "auc_results.xlsx")


 # --------------------------------------- Friedman rank sum test ----------------------------------

auc_results_matrix <- as.matrix(auc_results)
friedman.test(auc_results_matrix)



# --------------------------------------- Nemenyi test --------------------------------------------

###### The Nemenyi test is not working... ask for an alternative test to be used
auc <- posthoc.friedman.nemenyi.test(auc_results_matrix)
summary(auc)

## Alternate test
#library(scales)
#library(tibble)
#library(NSM3)
#pWNMT(auc_results_matrix, method = "Asymptotic")

#auc <- frdAllPairsNemenyiTest(auc_results_matrix)
#summary(auc)

auc$p.value[9, ]


ibs_results <- data.frame("cox" = res.cox$ibs, 
                          "aft" = res.aft$ibs,
                          "rsf" = res.rsf$ibs,
                          "mtlr" = res.mtlr$ibs,
                          "ssmtlr" = res.ssmtlr$ibs)


# --------------------------------------------------- Friedman rank sum test -------------------------------------------------

ibs_results_matrix <- as.matrix(ibs_results)
friedman.test(ibs_results_matrix)


# --------------------------------------------------- Nemenyi test ------------------------------------------------
ibs <- posthoc.friedman.nemenyi.test(ibs_results_matrix)

summary(ibs)
ibs$p.value[9, ]
