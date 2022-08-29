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
library("tcltk")
library("randomForestSRC")
library("ranger")

options(scipen=200)
options(rf.cores = -1)


# ******************************************** seer *************************************************
#####################################################################################################
########################################### load data ###############################################
#####################################################################################################
data(Melanoma)
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

folds <- createFolds(train$status, k = 5)
gbm_grid <- expand.grid(nsplit = c(5, 10, 15, 20, 25), 
                        ntree = c(100, 200, 300, 400, 500), 
                        mtry = c(2, 3, 4),
                        nodesize = c(20, 30, 40, 50, 60))


##################### No.1  RSF ###############################
# Function for CV
gbm.cv <- function(data, nsplit, ntree, mtry, nodesize, IBS = FALSE){
  aucs <- c() # Initiate
  ibss <- c()
  
  for (ifold in 1:length(folds)){
    itr <- unlist(folds[-ifold])
    Xtrain <- data[itr,]
    Xvalid <- data[-itr,]
    # Fit gbm
    model <- rfsrc(Surv(time, status) ~ age + sex + event + invasion + ici
                   + ulcer + epicel + thick,
                   splitrule = "logrankCR",
                   data = Xtrain, nsplit = nsplit, ntree = ntree, 
                   mtry = mtry, nodesize = nodesize,
                   ntime = c(12, 24, 36, 48, 60, 72, 84, 96), tree.err = FALSE)
    
    # Predict on validation set
    pred <- predict(model, newdata = Xvalid, proximity = FALSE)
    auc <- 1 - pred$err.rate[nrow(pred$err.rate), 1]
    aucs <- c(aucs, auc)
    
    ibs.v <- NA
    if(IBS){
      # Using pec for IBS estimation
      ibs.v <- pec(object = model,
                   formula = Surv(time, status) ~ age + sex + event + invasion + ici
                   + ulcer + epicel + thick,
                   cens.model = "marginal", cause = 1, 
                   times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                   data = Xvalid, verbose = F, maxtime = 200)
      ibss <- c(ibss, list(ibs.v))
    }
  }
  # Return the cindex vector
  if(IBS){
    return(list("aucs" = aucs, "ibss" = ibss))
  }
  else{
    return(aucs)
  }
}

# Loop over the rf_grid to screen best hyperparameters
gbm.cv.results <- data.frame(gbm_grid)
gbm.cv.results$auc <- 0


for(ind in 1:dim(gbm_grid)[1]){
  nsplit <- gbm_grid[ind,1]
  ntree <- gbm_grid[ind,2]
  mtry <- gbm_grid[ind,3]
  nodesize <- gbm_grid[ind,4]
  
  # Run fivefold CV
  gbm.rets <- gbm.cv(data = train, nsplit = nsplit, ntree = ntree, mtry = mtry, 
                     nodesize = nodesize, IBS = FALSE)
  # Save the mean to rf.cv.results
  gbm.cv.results[ind, 5] <- mean(gbm.rets[[1]])
  
  # Write report
  cat("Model :: ", ind, "/", dim(gbm_grid)[1], "nsplit = ", nsplit, "ntree = ",
      ntree, "mtry = ", mtry,
      "nodesize = ", nodesize,
      ":: auc = ", gbm.cv.results[ind, 5], "\n")
  save(gbm.cv.results, file="gbm.cv.results.once.Rdata")
}

save(gbm.cv.results, file="gbm.cv.results.once.Rdata")
load(file="gbm.cv.results.once.Rdata")

ind.best <- which.max(gbm.cv.results$auc)


# -------------------------------------------------------------------------------------------
############################################# No.1 RSF ######################################
# -------------------------------------------------------------------------------------------
###################################### modeling ##################################

model <- rfsrc(Surv(time, status) ~ age + sex + event + invasion + ici
               + ulcer + epicel + thick,
               splitrule="logrankCR",
               data = train, 
               nsplit = gbm.cv.results[ind.best, 1], 
               ntree = gbm.cv.results[ind.best, 2], 
               mtry = gbm.cv.results[ind.best, 3],
               nodesize = gbm.cv.results[ind.best, 4],
               ntime = c(12, 24, 36, 48, 60, 72, 84, 96),
               tree.err = TRUE)

save(model, file = "model_rsf.R")
load(file = "model_rsf.R")
print(model)


####################################### cindex ###################################
cindex.train.rsf <- 1 - model$err.rate[nrow(model$err.rate), 1]
cindex.train.rsf

pred.rsf <- predict(model, newdata = test, proximity = FALSE, outcome = "test")
pred.rsf

cindex.test.rsf <- 1 - pred.rsf$err.rate[nrow(pred.rsf$err.rate), 1]
cindex.test.rsf


####################################### ibs ###################################
Melanoma <- na.omit(Melanoma)

# Using pec for IBS estimation
ibs.train.rsf <- pec(object = model,
                     formula = Surv(time, status) ~ age + sex + event + invasion + ici
                     + ulcer + epicel + thick,
                     cens.model = "marginal", cause = 1, 
                     times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                     data = train, verbose = F, maxtime = 5000)
ibs.train.v.rsf <- crps(ibs.train.rsf)[2]
ibs.train.v.rsf


ibs.test.rsf <- pec(object = model,
                    formula = Surv(time, status) ~ age + sex + event + invasion + ici
                    + ulcer + epicel + thick,
                    cens.model = "marginal", cause = 1, 
                    times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                    data = test, verbose = F, maxtime = 5000)
ibs.test.v.rsf <- crps(ibs.test.rsf)[2]
ibs.test.v.rsf



####################################### bootstrap ###################################
Beta.rsf <- function(data, indices){
  d <- data[indices,]
  model.rsf <- rfsrc(Surv(time, status) ~ age + sex + event + invasion + ici
                     + ulcer + epicel + thick,
                     splitrule = "logrankCR", 
                     data = d, 
                     nsplit = gbm.cv.results[ind.best, 1], 
                     ntree = gbm.cv.results[ind.best, 2], 
                     mtry = gbm.cv.results[ind.best, 3],
                     nodesize = gbm.cv.results[ind.best, 4],
                     ntime = c(12, 24, 36, 48, 60, 72, 84, 96),
                     tree.err = TRUE)
  cindex.train <- 1 - model.rsf$err.rate[nrow(model.rsf$err.rate), 1]
  
  pred.rsf <- predict(model.rsf, newdata = test, proximity = FALSE, outcome = "test")
  cindex.test <- 1 - pred.rsf$err.rate[nrow(pred.rsf$err.rate), 1]
  
  
  ibs.train <- pec(object = model.rsf,
                   formula = Surv(time, status) ~ age + sex + event + invasion + ici
                   + ulcer + epicel + thick,
                   cens.model = "marginal", cause = 1, 
                   times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                   data = d, 
                   verbose = F, maxtime = 1000)
  ibs.train.v <- crps(ibs.train)[2]
  
  
  ibs.test <- pec(object = model.rsf,
                  formula = Surv(time, status) ~ age + sex + event + invasion + ici
                  + ulcer + epicel + thick,
                  cens.model = "marginal", cause = 1,
                  times = c(12,24,36,48,60,72,84,96), exact = TRUE,
                  data = test, 
                  verbose = F, maxtime = 5000)
  ibs.test.v <- crps(ibs.test)[2]
  return(c(cindex.train, cindex.test, ibs.train.v, ibs.test.v))
}

result.ci.rsf <- boot(data = train, statistic = Beta.rsf, R = 100)
save(result.ci.rsf, file = "result.ci.rsf.R")
load(file = "result.ci.rsf.R")


round(mean(result.ci.rsf$t[, 1]), 4)
round(mean(result.ci.rsf$t[, 2]), 4)
round(mean(result.ci.rsf$t[, 3]), 4)
round(mean(result.ci.rsf$t[, 4]), 4)


get_cis <- function(idx){
  boots.train.cindex.rsf <- result.ci.rsf
  boots.train.cindex.rsf$t0 <- result.ci.rsf$t0[idx]
  boots.train.cindex.rsf$t <- as.matrix(result.ci.rsf$t[, idx])
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
