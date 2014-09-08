#-------------------------------------------------------------------------------
# This R script is to perform model selections with and without offset using
# 1) GLM Poisson, 2) GLM Negative Binomial, 3) Rpart Poisson, 4) glmnet Poisson
# 5) gbm POisson. Our goal is to identify important variables on modeling 
# accident frequency using different models. Step() and variable importance will
# be used.
#-------------------------------------------------------------------------------

library(MASS)
library(rpart)
library(randomForest)
library(glmnet)
library(gbm)
library(texreg)
library(ggplot2)

#-- Load data
load("data/accidents.RData")

#-- Define formulam. Put back Lat and Long
fmla <- as.formula("X5YrCrashCount ~ AreaType +
                            IntCat + IntTCType + LegRtType + LegSpeed + LegTCType + LegType + 
                            LegWidth + Lighting + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                            MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                            NumLanes + NumSegs + Offset + OffsetDist + 
                            OneWay + PaveType + PedCross + RTChannel + 
                            RTLanes + RTLnLength + RTMoveCtrl + RTWidth + 
                            Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                            TotalAADT + log(TotalT5YearInM) + TurnProhib + Lat + Long") # keep log(TotalT5YearInM) here as a predictor for regression model

fmla.offset <- as.formula("X5YrCrashCount ~ AreaType +
                             IntCat + IntTCType + LegRtType + LegSpeed + LegTCType + LegType + 
                             LegWidth + Lighting + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                             MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                             NumLanes + NumSegs + Offset + OffsetDist + 
                             OneWay + PaveType + PedCross + RTChannel + 
                             RTLanes + RTLnLength + RTMoveCtrl + RTWidth + 
                             Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                             TotalAADT + TurnProhib + offset(log(TotalT5YearInM)) + Lat + Long")

#-------------------------------------------------------------------------------
# Candidate models are built in this section
#-------------------------------------------------------------------------------

#-------------#
#-- Poisson --#
#-------------#
poisReg <- glm(fmla, data = a, family = "poisson")
poisReg.offset <- glm(fmla.offset, data = a, family = "poisson")

# step.poisReg <- step(poisReg)
# step.poisReg.offset <- step(poisReg.offset)

#-- Quisi Poisson (No stepwise since AIC is not defined) --#
quasiPois <- glm(fmla, data = a, family = "quasipoisson")
quasiPois.offset <- glm(fmla.offset, data = a, family = "quasipoisson")

#-----------------------#
#-- Negative Binomial --#
#-----------------------#
negBino <- glm.nb(fmla, data = a)
negBino.offset <- glm.nb(fmla.offset, data = a)

# step.negBino <- step(negBino)
# step.negBino.offset <- step(negBino.offset)

#-------------------#
#-- rpart Poisson --#
#-------------------#
treePois <- rpart(fmla, data = a, method = "poisson")
treePois.offset <- rpart(fmla.offset, data = a, method = "poisson")

par(oma=c(2,3,2,2))
barplot(treePois$variable.importance, horiz = T, las = 1)
title("Variable Importance Poisson Regression Trees without Offset")

barplot(treePois.offset$variable.importance, horiz = T, las = 1)
title("Variable Importance Poisson Regression Trees with Offset")

#--------------------#
#-- glmnet poisson --#
#--------------------#
#- process the x, y for glmnet. y = X5YrCrashCount, x is from the flma.string
fmla.string <- "AreaType + IntCat + IntTCType + LegRtType + LegSpeed + LegTCType + LegType + LegWidth + Lighting + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + MedWidth + MergeLanes + NextPIDist + NumberLegs + NumLanes + NumSegs + Offset + OffsetDist + OneWay + PaveType + PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + TotalAADT + TurnProhib + Lat + Long" # put log(TotalT5YearInM) in the model matrix
fmla.string <- strsplit(fmla.string, " \\+ ")[[1]]

#- glmnet poisson no offset with LASSO
#- categorical variables need to be treated http://statweb.stanford.edu/~tibs/lasso/fulltext.pdf
# glmnetPois.lasso <- glmnet(model.matrix(fmla, data = a), as.matrix(a[, "X5YrCrashCount"]), family = "poisson") # X5YrCrashCount is in the model matrix?
# glmnetPois.lasso2 <- glmnet(as.matrix(model.frame(fmla, data = a)), as.matrix(a[, "X5YrCrashCount"]), family = "poisson")

# Let's stick to using original data first. Treat them as continuous variables
glmnetPois.lasso <- glmnet(as.matrix(cbind(a[, fmla.string], log(a[, "TotalT5YearInM"]))), as.matrix(a[, "X5YrCrashCount"]), 
                                  family = "poisson", offset = as.matrix(log(a[, "TotalT5YearInM"])))

par(oma=c(2,3,2,4))
plot(glmnetPois.lasso)
title("glmnet Poisson LASSO without Offset", line = 3)
vnat <- coef(glmnetPois.lasso)
vnat <- vnat[-1,ncol(vnat)] # remove the intercept term
axis(4, at = vnat, line = -0.5, label = c(fmla.string, "log(TotalT5YearInM)"), las = 1, tick = FALSE, cex.axis = 0.7) # problem with categorial dummies

#- glmnet poisson offset with LASSO
glmnetPois.offset.lasso <- glmnet(as.matrix(a[, fmla.string]),as.matrix(a[, "X5YrCrashCount"]), 
                                  family = "poisson", offset = as.matrix(log(a[, "TotalT5YearInM"])))

plot(glmnetPois.offset.lasso)
title("glmnet Poisson LASSO with Offset", line = 3)
vnat <- coef(glmnetPois.offset.lasso)
vnat <- vnat[-1,ncol(vnat)] # remove the intercept term
axis(4, at = vnat, line = -0.5, label = fmla.string, las = 1, tick = FALSE, cex.axis = 0.7)

#--------------------#
#-- gbm poisson --#
#--------------------#
gbmPois <- gbm(fmla, data = a, distribution = "poisson")
summary(gbmPois, main = "ssd")
gbmPois.offset <- gbm(fmla.offset, data = a, distribution = "poisson")
summary(gbmPois.offset)


#-------------------------------------------------------------------------------
# In this section, we summarized useful variables in each model (pois, nb)
#-------------------------------------------------------------------------------
alpha <- 0.05
summary(poisReg)$coefficients[summary(poisReg)$coefficients[, 4] < alpha, ]
summary(poisReg.offset)$coefficients[summary(poisReg.offset)$coefficients[, 4] < alpha, ]
summary(negBino)$coefficients[summary(negBino)$coefficients[, 4] < alpha, ]
summary(negBino.offset)$coefficients[summary(negBino.offset)$coefficients[, 4] < alpha, ]

barplot(sort(summary(poisReg)$coefficients[summary(poisReg)$coefficients[, 4] < alpha, 4], 
             decreasing = T), horiz = T, las = 1, main = "Poisson Regression no Offset",
        xlab = "p-value")
barplot(sort(summary(poisReg.offset)$coefficients[summary(poisReg.offset)$coefficients[, 4] < alpha, 4], 
             decreasing = T), horiz = T, las = 1, main = "Poisson Regression with Offset",
        xlab = "p-value")
barplot(sort(summary(negBino)$coefficients[summary(negBino)$coefficients[, 4] < alpha, 4], 
             decreasing = T), horiz = T, las = 1, main = "Negative Binomial no Offset",
        xlab = "p-value")
barplot(sort(summary(negBino.offset)$coefficients[summary(negBino.offset)$coefficients[, 4] < alpha, 4], 
             decreasing = T), horiz = T, las = 1, main = "Negative Binomial with Offset",
        xlab = "p-value")


#-------------------------------------------------------------------------------
# In this section, we scale p-values and variable importance for comparison
#-------------------------------------------------------------------------------
# 1) sort variable names
# 2) keep only one for categorical variables
# 3) Fill in zeros when using importances
allVariables <- c(fmla.string, "log(TotalT5YearInM)") # complete set of variables names
allVariables <- allVariables[order(allVariables)]
varImpStandard <- function(v) { # scale variable "importance" from 0 to 100
  v <- abs(v) * 100 / max(abs(v)) # use absolute values when only coeffcients are provided
  if (length(v) == length(allVariables)) return(v)
  else if (length(v) < length(allVariables)) { # deal with tree importance and non-offset
    s <- rep(0, length(allVariables))
    names(s) <- allVariables
    matchedVariables <- allVariables[which(allVariables %in% names(v))]
    s[matchedVariables] <- v[matchedVariables]
    return(s)  
  }
  else { # deal with dummy variables, exclude intercept
    f.temp <- function(x) { # strip the last digit for categorical variables
      if (substr(x, nchar(x), nchar(x)) %in% c(2:9)) return(substr(x, 1, nchar(x) - 1))
      else return(x)
    }
    v.names.stripped <- unlist(lapply(names(v), FUN = f.temp))
    v.agg <- aggregate(v, by = list(v.names.stripped), FUN = min)
    v <- v.agg$x
    names(v) <- v.agg$Group.1
    return(v)
  }    
}

#dump common 0 rows
varImp_1 <- varImpStandard(summary(poisReg)$coefficients[-1, 4])
varImp_2 <- varImpStandard(summary(poisReg.offset)$coefficients[-1, 4])
varImp_3 <- varImpStandard(summary(negBino)$coefficients[-1, 4])
varImp_4 <- varImpStandard(summary(negBino.offset)$coefficients[-1, 4])
varImp_5 <- varImpStandard(treePois$variable.importance)
varImp_6 <- varImpStandard(treePois.offset$variable.importance)
varImp_7 <- varImpStandard(glmnetPois.lasso$beta[, glmnetPois.lasso$dim[2]])
varImp_8 <- varImpStandard(glmnetPois.offset.lasso$beta[, glmnetPois.offset.lasso$dim[2]])

importances <- rbind(cbind(allVariables, varImp_1, "Poisson Regression"),
                     cbind(allVariables, varImp_2, "Offset Poisson Regression"),
                     cbind(allVariables, varImp_3, "Negative Binomia"),
                     cbind(allVariables, varImp_4, "Offset Negative Binomial"),
                     cbind(allVariables, varImp_5, "Regression Tree"),
                     cbind(allVariables, varImp_6, "Offset Regression Tree"),
                     cbind(allVariables, varImp_7, "glmnet Poisson"),
                     cbind(allVariables, varImp_8, "Offset glmnet Poisson"))                                              )

importances <- cbind(varImpStandard(summary(poisReg)$coefficients[, 4]),
                     varImpStandard(summary(poisReg.offset)$coefficients[, 4]),
                     varImpStandard(summary(negBino)$coefficients[, 4]),
                     varImpStandard(summary(negBino.offset)$coefficients[, 4]),
                     varImpStandard(treePois$variable.importance),
                     varImpStandard(treePois.offset$variable.importance),
                     varImpStandard(glmnetPois.lasso$beta[, glmnetPois.lasso$dim[2]]),
                     varImpStandard(glmnetPois.offset.lasso$beta[, glmnetPois.offset.lasso$dim[2]])
                     )
names(importances) <- c("Poisson Regression", "Offset Poisson Regression", 
                        "Negative Binomial", "Offset Negative Binomial",
                        "Regression Tree", "Offset Regression Tree",
                        "glmnet Poisson", "Offset glmnet Poisson")

barplot(varImpStandard(s), horiz = T, las = 1, main = "Poisson Regression no Offset",
        xlab = "p-value")


#-------------------------------------------------------------------------------
# In this section, we tend to evaluate the models we obtained
#-------------------------------------------------------------------------------



