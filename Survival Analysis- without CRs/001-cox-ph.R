setwd("..../Survival analysis without CR/")
rm(list=ls())

library(survival)
library(foreign)
library(rms)
library(ROCR)
library(caret)
library(cluster)
library(randomForestSRC)
library(parallel)
library(cmprsk)
library(aod)
library(pec)
library(gbm)
library(ranger)
library(boot)
options(scipen=200)
options(contrasts=c("contr.treatment", "contr.treatment"))


# ******************************************** seer *************************************************
#####################################################################################################
########################################### load data ###############################################
#####################################################################################################
data <- read.csv("Birth Control.csv")
str(data)
#names(data)

sum(is.na(data))
#cancer <- na.omit(cancer)
#sum(is.na(cancer))

# renamng the column
colnames(data)[colnames(data) == "status"] <- "os"

## pre-processing
#data$os <- ifelse(data$os == 2, 1, 0)
data$homeStyle <- ifelse(data$homeStyle == "urban", 1, 0)

str(data)

# Splitting the dataset
library(caTools)
split <- sample.split(data$time, SplitRatio = 0.8)

train <- subset(data, split == TRUE)
test <- subset(data, split == FALSE)



# -------------------------------------------------------------------------------------------
##################################### No.1 Cox survival #####################################
# -------------------------------------------------------------------------------------------

###################################### modeling ##################################
model <- coxph(Surv(time, os) ~ age + region + homeStyle + wealth + children,
               data = train, x = TRUE)


####################################### cindex ###################################
concordance(model)
concordance(model, newdata = test)



####################################### ibs ###################################
# Using pec for IBS estimation
ibs.train <- pec(object = model,
                 formula = Surv(time, os) ~ age + region + homeStyle + wealth + children,
                 cens.model = "marginal", times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                 data = train, verbose = F, maxtime = 1000)
ibs.train.v <- crps(ibs.train)[2]
ibs.train.v


ibs.test <- pec(object = model,
                formula = Surv(time, os) ~ age + region + homeStyle + wealth + children,
                cens.model = "marginal", times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                data = test, verbose = F, maxtime = 1000)
ibs.test.v <- crps(ibs.test)[2]
ibs.test.v


####################################### bootstrap ###################################
Beta.model <- function(data, indices){
  d <- data[indices,]
  model <- coxph(Surv(time, os) ~ age + region + homeStyle + wealth + children,
                 data = d,
                 x = TRUE)
  cindex.train <- concordance(model)$concordance
  
  cindex.test <- concordance(model, newdata = test)$concordance
  
  
  ibs.train <- pec(object = model,
                   formula = Surv(time, os) ~ age + region + homeStyle + wealth + children,
                   cens.model = "marginal", times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                   data = d,
                   verbose = F, maxtime = 1000)
  ibs.train.v <- crps(ibs.train)[2]
  
  
  ibs.test <- pec(object = model,
                  formula = Surv(time, os) ~ age + region + homeStyle + wealth + children,
                  cens.model = "marginal", times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                  data = test, verbose = F, maxtime = 1000)
  ibs.test.v <- crps(ibs.test)[2]
  
  return(c(cindex.train, cindex.test, ibs.train.v, ibs.test.v))
}

result.ci <- boot(data = train, statistic = Beta.model, R = 100)
save(result.ci, file = "result.ci.cox.R")
load(file = "result.ci.cox.R")


round(mean(result.ci$t[, 1]), 4)
round(mean(result.ci$t[, 2]), 4)
round(mean(result.ci$t[, 3]), 4)
round(mean(result.ci$t[, 4]), 4)


get_cis <- function(idx){
  boots.train.cindex.rsf <- result.ci
  boots.train.cindex.rsf$t0 <- result.ci$t0[idx]
  boots.train.cindex.rsf$t <- as.matrix(result.ci$t[, idx])
  res <- boot.ci(boots.train.cindex.rsf, conf = 0.95, type = "perc")
  return(res$percent[4:5])
}
# train cindex
round(get_cis(1), 4)

# test cindex
round(get_cis(2), 4)


# train ibs
round(get_cis(3), 4)

# test ibs
round(get_cis(4), 4)

