library(MASS)
library(rpart)
library(randomForest)

#-- load data
rm(list = ls())
load("data/accidents.RData") # dataframe is a

#-- define two formula, removed variables X5YrCrashCount, IntID, LegID, Lat, Long, X1YrCrashCount, City, County, TotalAADT
#-- Remove all ADDT variables? either in offset or appear once as predictor?
fmla.OffsetNo <- as.formula("X5YrCrashCount ~ AreaType + IntAADT +
                            IntCat + IntTCType + LegRtType + LegSpeed + LegTCType + LegType + 
                            LegWidth + Lighting + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                            MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                            NumLanes + NumSegs + Offset + OffsetDist + 
                            OneWay + PaveType + PedCross + RTChannel + 
                            RTLanes + RTLnLength + RTMoveCtrl + RTWidth + 
                            Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                            TotalAADT + TotalT5YearInM + TurnProhib")

fmla.OffsetYes <- as.formula("X5YrCrashCount ~ AreaType + IntAADT +
                             IntCat + IntTCType + LegRtType + LegSpeed + LegTCType + LegType + 
                             LegWidth + Lighting + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                             MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                             NumLanes + NumSegs + Offset + OffsetDist + 
                             OneWay + PaveType + PedCross + RTChannel + 
                             RTLanes + RTLnLength + RTMoveCtrl + RTWidth + 
                             Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                             TotalAADT + TurnProhib + offset(log(TotalT5YearInM))") # offset(log(TotalT5YearInM))


#-- Poisson Regression
# diagnostics of possion no offset
ps.OffSetNo <- glm(fmla.OffsetNo, data = a, family = "poisson")

sum(residuals(ps.OffSetNo, type = "pearson")^2) # pearson residuals for goodness of fit

(1 - exp((ps.OffSetNo$dev - ps.OffSetNo$null) / dim(a)[1])) / (1 - exp(-ps.OffSetNo$null / dim(a)[1])) # Naglekerke R-square

obj <- glmulti(fmla.OffsetNo, data = a, family = "poisson", level = 1, crit = bic) # use glmulti package for glm model selection

step(ps.OffSetNo) # step-wise model selection

ps.OffSetYes <- glm(fmla.OffsetYes, data = a, family = "poisson")

sum(residuals(ps.OffSetYes, type = "pearson")^2) # pearson residuals for goodness of fit

(1 - exp((ps.OffSetYes$dev - ps.OffSetYes$null) / dim(a)[1])) / (1 - exp(-ps.OffSetYes$null / dim(a)[1])) # Naglekerke R-square

glmulti(fmla.OffsetYes, data = a, family = "poisson", level = 1, maxit = 20) # use glmulti package for glm model selection

step(ps.OffSetYes) # step-wise model selection


#-- Quisi Poisson
qp.OffSetNo <- glm(fmla.OffsetNo, data = a, family = "quasipoisson") # no complaints using quasi-poisson
qp.OffSet <- glm(fmla.OffsetYes, data = a, family = "quasipoisson")

#-- overdispersed, negative binomial
nb.OffSetNo <- glm.nb(fmla.OffsetNo, data = a)
nb.OffSetYes <- glm.nb(fmla.OffsetYes, data = a)

#-- Hurdle Model
#-- Zero-inated models
#-- Regression Tree. Do this work as offset? Same result
rt.OffSetNo <- rpart(fmla.OffsetNo, data = a, method = "poisson")
rt.OffSetYes <- rpart(fmla.OffsetYes, data = a, method = 'poisson')

#-- Random Forests
rf.OffSetNo <- randomForest(fmla.OffsetNo, data = a)
varImpPlot(rf.OffSetNo) # variable importance plot

rf.OffSetYes <- randomForest(fmla.OffsetYes, data = a) # somehow it works mystery
varImpPlot(rf.OffSetYes) # variable importance plot






