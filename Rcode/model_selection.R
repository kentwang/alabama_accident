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
# library(texreg)
library(ggplot2)

#-- Load data
load("data/accidents.RData")

#-- Define formulam. Put back Lat and Long
fmla <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
                    LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
                    LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                    MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                    NumLanes + Offset + OffsetDist + OneWay + PaveType + 
                    PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
                    RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                    log(Traffic) + TurnProhib") 

fmla.offset <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
                    LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
                    LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                    MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                    NumLanes + Offset + OffsetDist + OneWay + PaveType + 
                    PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
                    RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                    offset(log(Traffic)) + TurnProhib") # 38 including the offset


#-------------------------------------------------------------------------------
# Candidate models are built in this section
#-------------------------------------------------------------------------------

#-------------#
#-- Poisson --#
#-------------#
poisReg <- step(glm(fmla, data = a, family = "poisson"))
poisReg.offset <- step(glm(fmla.offset, data = a, family = "poisson"))

# step.poisReg <- step(poisReg)
# step.poisReg.offset <- step(poisReg.offset)

#-- Quisi Poisson (No stepwise since AIC is not defined) --#
quasiPois <- glm(fmla, data = a, family = "quasipoisson")
quasiPois.offset <- glm(fmla.offset, data = a, family = "quasipoisson")

#-----------------------#
#-- Negative Binomial --#
#-----------------------#
negBino <- step(glm.nb(fmla, data = a))
negBino.offset <- step(glm.nb(fmla.offset, data = a))

# step.negBino <- step(negBino)
# step.negBino.offset <- step(negBino.offset)

#-------------------#
#-- rpart Poisson --#
#-------------------#
treePois <- rpart(fmla, data = a, method = "poisson")
treePois.offset <- rpart(fmla.offset, data = a, method = "poisson")

op <- par(mar = c(1, 0, 1, 0))
plot(treePois, main = "Poisson Regression Tree no Offset")
text(treePois, use.n=TRUE, cex = .8)
par(op)


#--------------------#
#-- glmnet poisson --#
#--------------------#

#- glmnet poisson no offset with LASSO
#- categorical variables need to be treated http://statweb.stanford.edu/~tibs/lasso/fulltext.pdf
# glmnetPois.lasso2 <- glmnet(as.matrix(model.frame(fmla, data = a)[,-1]), as.matrix(a[, "X5YrCrashCount"]), family = "poisson")

glmnetPois.cv <- cv.glmnet(model.matrix(fmla, data = a)[,-1], as.matrix(a[, "X5YrCrashCount"]), family = "poisson") # X5YrCrashCount is in the model matrix?
glmnetPois.cv.offset <- cv.glmnet(model.matrix(fmla.offset, data = a)[,-1], as.matrix(a[, "X5YrCrashCount"]), 
                               family = "poisson", offset = log(a$Traffic))

std = apply(model.matrix(fmla, data = a)[,-1],2,sd)
beta <- as.matrix(coef(glmnetPois.cv, s="lambda.1se"))
beta.names <- rownames(beta)
beta <- as.numeric(beta)
names(beta) <- beta.names

# Todo: Some issue with selecting best lambda for offset cv.glmnet. Manually select
beta.offset <- as.matrix(coef(glmnetPois.cv.offset, s="lambda.1se"))
beta.offset.names <- rownames(beta.offset)
beta.offset <- as.numeric(beta.offset)
names(beta.offset) <- beta.offset.names
## importance by abs(beta); pick from factor


# Let's stick to using original data first. Treat them as continuous variables
glmnetPois.lasso <- glmnet(model.matrix(fmla, data = a)[,-1], as.matrix(a[, "X5YrCrashCount"]), family = "poisson")

par(oma=c(2,3,2,4))
plot(glmnetPois.lasso)
title("glmnet Poisson LASSO without Offset", line = 3)
vnat <- coef(glmnetPois.lasso)
vnat <- vnat[-1,ncol(vnat)] # remove the intercept term
axis(4, at = vnat, line = -0.5, label = rownames(glmnetPois.lasso$beta), las = 1, tick = FALSE, cex.axis = 0.7) # problem with categorial dummies

abline(v = sum(abs(beta)), col = "blue", lwd = 2, lty = "dotted")

#- glmnet poisson offset with LASSO
glmnetPois.offset.lasso <- glmnet(model.matrix(fmla.offset, data = a)[,-1], as.matrix(a[, "X5YrCrashCount"]), 
                                  family = "poisson", offset = log(a$Traffic))

par(oma=c(1,1,1,4))
plot(glmnetPois.offset.lasso)
title("glmnet Poisson LASSO with Offset", line = 3)
vnat <- coef(glmnetPois.offset.lasso)
vnat <- vnat[-1,ncol(vnat)] # remove the intercept term
axis(4, at = vnat, line = -0.5, label = rownames(glmnetPois.offset.lasso$beta), las = 1, tick = FALSE, cex.axis = 0.7)

abline(v = sum(abs(beta.offset)), col = "blue", lwd = 2, lty = "dotted")
text(x = 12, y = -6, "CV Optimal Coef")
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
fmla.string <- strsplit(as.character(fmla)[3], " \\+ ")[[1]]
# allVariables <- c(fmla.string, "log(TotalT5YearInM)") # complete set of variables names
allVariables <- fmla.string[order(fmla.string)]


#dump common 0 rows, but there are none
varImp_1 <- varImpStandard(summary(poisReg)$coefficients[-1, 3])
varImp_2 <- varImpStandard(summary(poisReg.offset)$coefficients[-1, 3])
varImp_3 <- varImpStandard(summary(negBino)$coefficients[-1, 3])
varImp_4 <- varImpStandard(summary(negBino.offset)$coefficients[-1, 3])
varImp_5 <- varImpStandard(treePois$variable.importance)
varImp_6 <- varImpStandard(treePois.offset$variable.importance)
varImp_7 <- varImpStandard(beta[-1])
varImp_8 <- varImpStandard(beta.offset[-1]) #glmnetPois.offset.lasso selected nothing?
varImp <- as.vector(c(varImp_1, varImp_2, varImp_3, varImp_4, 
                          varImp_5, varImp_6, varImp_7, varImp_8))

model.names <- c("Poisson Regression", "Poisson Regression Offset",
                 "Negative Binomial","Negative Binomial Offset",
                 "Regression Tree", "Regression Tree Offset",
                 "glmnet Poisson", "glmnet Poisson Offset")

importances <- as.data.frame(cbind(rep(allVariables, length(model.names)),
                                   rep(model.names, each = length(allVariables))))
importances$importance <- varImp
names(importances) <- c("var", "model", "importance")

# All together is kinda fine
ggplot(importances, aes(factor(var), importance, fill = model)) + 
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip()

# offset
ggplot(importances[grep("Offset", importances$model), ], aes(factor(var), importance, fill = model)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(text = element_text(size = 16)) + 
  coord_flip()

# no offset
ggplot(importances[-grep("Offset", importances$model), ], aes(factor(var), importance, fill = model)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(text = element_text(size = 16)) + 
  coord_flip()

#-------------------------------------------------------------------------------
# In this section, we tend to evaluate the models we obtained
#-------------------------------------------------------------------------------


