setwd("C:/Users/USR/Documents/External Projects/MSc Survival analysis with deep learning/Survival analysis with CR/")

getwd()

rm(list=ls())




library("survival")
library("foreign")
library("rms")
library("ROCR")
library("DMwR2")
library("cluster")
library("riskRegression")
library("caret")
library("parallel")
library("cmprsk")
library("aod")
library("pec")
library("gbm")
library("boot")
library("prodlim")


options(scipen=200)
options(rf.cores = -1)


# ******************************************** seer *************************************************
#####################################################################################################
########################################### load data ###############################################
#####################################################################################################


data("Melanoma")
str(Melanoma)
#attach(Melanoma)
names(Melanoma)

sum(is.na(Melanoma))
#Melanoma <- na.omit(Melanoma)
#sum(is.na(Melanoma))



# Splitting the dataset
library(caTools)
split <- sample.split(Melanoma$time, SplitRatio = 0.8)

train <- subset(Melanoma, split == TRUE)
test <- subset(Melanoma, split == FALSE)
str(train)

################################### No.1  Fine-Gray competing risks #################################

model <- FGR(Hist(time, status) ~ age + sex + event + invasion + ici
             + ulcer + epicel + thick,
             data = train,
             cause = 1)
summary(model)

cindex.train  <- pec::cindex(model,
                             formula = Hist(time, status) ~ age + sex + event + invasion + ici
                             + ulcer + epicel + thick,
                             data = train,
                             cens.model = "marginal", 
                             cause = 1, eval.times = 1000)
cindex.train <- cindex.train$AppCindex$FGR
cindex.train

cindex.test  <- pec::cindex(model,
                            formula = Hist(time, status) ~ age + sex + event + invasion + ici
                            + ulcer + epicel + thick,
                            data = test,
                            cens.model = "marginal", cause = 1,
                            eval.times = 1000)
cindex.test <- cindex.test$AppCindex$FGR
cindex.test



ibs.train <- pec(object = model,
                 formula = Hist(time, status) ~ age + sex + event + invasion + ici
                 + ulcer + epicel + thick,
                 cens.model = "marginal", cause = 1, exact = TRUE,
                 data = train, verbose = F)
ibs.v.train <- crps(ibs.train)[2]
ibs.v.train

ibs.test <- pec(object = model,
                formula = Hist(time, status) ~ age + sex + event + invasion + ici
                + ulcer + epicel + thick,
                cens.model = "marginal", cause = 1, exact = TRUE,
                data = test, verbose = F)
ibs.v.test <- crps(ibs.test)[2]
ibs.v.test


####################################### bootstrap ###################################
Beta.model <- function(data, indices){
  d <- data[indices,]
  model <- FGR(Hist(time, status) ~ age + sex + event + invasion + ici
               + ulcer + epicel + thick,
               data = d,
               cause = 1)
  cindex.train  <- pec::cindex(model,
                               formula = Hist(time, status) ~ age + sex + event + invasion + ici
                               + ulcer + epicel + thick,
                               data = d,
                               cens.model = "marginal", cause = 1,
                               eval.times = 1000)
  
  cindex.train <- cindex.train$AppCindex$FGR
  
  cindex.test  <- pec::cindex(model,
                              formula = Hist(time, status) ~ age + sex + event + invasion + ici
                              + ulcer + epicel + thick,
                              data = test,
                              cens.model = "marginal", cause = 1,
                              eval.times = 1000)
  
  cindex.test <- cindex.test$AppCindex$FGR
  
  
  ibs.train <- pec(object = model,
                   formula = Hist(time, status) ~ age + sex + event + invasion + ici
                   + ulcer + epicel + thick,
                   cens.model = "marginal", cause = 1, exact = TRUE,
                   data = d,
                   verbose = F)
  ibs.v.train <- crps(ibs.train)[2]
  
  ibs.test <- pec(object = model,
                  formula = Hist(time, status) ~ age + sex + event + invasion + ici
                  + ulcer + epicel + thick,
                  cens.model = "marginal", cause = 1, exact = TRUE,
                  data = test, 
                  verbose = F)
  ibs.v.test <- crps(ibs.test)[2]
  return(c(cindex.train, cindex.test, ibs.v.train, ibs.v.test))
}

result.ci <- boot(data = train, statistic = Beta.model, R = 100)
result.ci

save(result.ci, file = "result.ci.fgr.R")
load(file = "result.ci.fgr.R")

mean(result.ci$t[,1])
round(mean(result.ci$t[,1]), 4)
round(mean(result.ci$t[,2]), 4)
round(mean(result.ci$t[,3]), 4)
round(mean(result.ci$t[,4]), 4)


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

#cancer <- data.frame(cancer)
#cancer

#library(writexl)
#write_xlsx(cancer, "Cancer.xlsx")
