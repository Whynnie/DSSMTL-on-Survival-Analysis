setwd("...../Survival analysis without CR/")

getwd()
rm(list = ls())

library(survival)
data = read.csv("..../Colon cancer.csv", header=TRUE)
str(data)
attach(data)
sex = as.factor(sex)
obstruct = as.factor(obstruct)
perfor = as.factor(perfor)
adhere = as.factor(adhere)
differ = ifelse(is.na(differ),1,differ)
differ = as.factor(differ)
extent = as.factor(extent)
surg=as.factor(surg)
node4 = as.factor(node4)
data = data.frame(time,	status, age,	sex,	obstruct,	perfor, adhere, node4,
                  differ, extent, surg)


##########  Parametrics Method #######################################
######################################################################

library(randomForest)
library(randomForestSRC)

library(pec)
library(rms)


## Variable Important
library(ggRandomForests)
Randomsurv <- rfsrc(Surv(time, status)~., data = data, nsplit = 1000, 
                    na.action = "na.impute")
varsel_pbc <- var.select(Randomsurv)

# ---------------------------------------- Without CRs --------------------------------------
plot(gg_vimp(Randomsurv),main="Variable Important")


#----------------------------------------- With CRs --------------------------------
#plot(gg_vimp(vimp.rfsrc(Randomsurv)),main="Variable Important")
