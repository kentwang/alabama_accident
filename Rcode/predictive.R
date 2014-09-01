#####################################################################
#  This script is to build predictive models on count of accidents.
#  Models to be used are Poisson Regression, Negative Binomial, Quasi-Poisson
#  regression tree in 'rpart' and Random Forest
#
#  Author: Ketong Wang
#####################################################################

library(MASS)
library(rpart)
library(randomForest)

#-- load data
rm(list = ls())
load("data/accidents.RData") # dataframe is a

#-- All variable names
# paste(sort(names(a)), collapse = " + ")
# "AADT + AreaType + City + County + IntAADT + IntCat + IntID + IntTCType + Lat + LegID + LegRtType + LegSpeed + LegTCType + LegType + LegWidth + Lighting + Long + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + MedWidth + MergeLanes + NextPIDist + NumberLegs + NumLanes + NumSegs + Offset + OffsetDist + OneWay + PaveType + PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + TotalAADT + TotalT5YearInM + TurnProhib + X1YrCrashCount + X5YrCrashCount"

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
ps.OffSetNo <- glm(fmla.OffsetNo, data = a, family = "poisson")
plot(predict(ps.OffSetNo), residuals(ps.OffSetNo), xlab = "Fitted", ylab = "Residuals")


ps.OffSetYes <- glm(fmla.OffsetYes, data = a, family = "poisson")

#-- Quisi Poisson
qp.OffSetNo <- glm(fmla.OffsetNo, data = a, family = "quasipoisson") # no complaints using quasi-poisson
# sum(residuals(ps.OffSetNo, type = "pearson")^2) # pearson residuals for goodness of fit
qp.OffSet <- glm(fmla.OffsetYes, data = a, family = "quasipoisson")

#-- overdispersed, negative binomial
bn.OffSetNo <- glm.nb(fmla.OffsetNo, data = a)
bn.OffSetYes <- glm.nb(fmla.OffsetYes, data = a)

#-- Hurdle Model

#-- Zero-inated models

#-- Regression Tree. Do this work as offset?
rt.OffSetNo <- rpart(fmla.OffsetNo, data = a, method = "poisson")
rt.OffSetYes <- rpart(fmla.OffsetYes, data = a, method = 'poisson')
# rt.OffSetYes <- rpart("cbind(TotalT5YearInM, X5YrCrashCount) ~ AreaType + IntAADT +
#                             IntCat + IntTCType + LegRtType + LegSpeed + LegTCType + LegType + 
#                             LegWidth + Lighting + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
#                             MedWidth + MergeLanes + NextPIDist + NumberLegs + 
#                             NumLanes + NumSegs + Offset + OffsetDist + 
#                             OneWay + PaveType + PedCross + RTChannel + 
#                             RTLanes + RTLnLength + RTMoveCtrl + RTWidth + 
#                             Rumble + SightLt + SightRt + SkewAngle + Terrain + 
#                             TotalAADT + TurnProhib", data = a, method = "poisson")

#-- Random Forests
rf.OffSetNo <- randomForest(fmla.OffsetNo, data = a)
varImpPlot(rf.OffSetNo) # variable importance plot

rf.OffSetYes <- randomForest(fmla.OffsetYes, data = a) # somehow it works mystery
varImpPlot(rf.OffSetYes) # variable importance plot







