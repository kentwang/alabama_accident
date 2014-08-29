#####################################################################
#  This script is to build predictive models on count of accidents.
#  Models to be used are Poisson Regression, Negative Binomial, regression
#  tree in 'rpart' and Random Forest
#
#  Author: Ketong Wang
#####################################################################

#-- load data
rm(list = ls())
load("data/accidents.RData") # dataframe is a

#-- All variable names
# paste(sort(names(a)), collapse = " + ")
# "AADT + AreaType + City + County + IntAADT + IntCat + IntID + IntTCType + Lat + LegID + LegRtType + LegSpeed + LegTCType + LegType + LegWidth + Lighting + Long + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + MedWidth + MergeLanes + NextPIDist + NumberLegs + NumLanes + NumSegs + Offset + OffsetDist + OneWay + PaveType + PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + TotalAADT + TotalT5YearInM + TurnProhib + X1YrCrashCount + X5YrCrashCount"

#-- define two formula, removed variables X5YrCrashCount, IntID, LegID, Lat, Long, X1YrCrashCount
fmla.OffsetNo <- as.formula("X5YrCrashCount ~ AADT + AreaType + City + County + IntAADT +
                            IntCat + IntTCType + LegRtType + LegSpeed + LegTCType + LegType + 
                            LegWidth + Lighting + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                            MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                            NumLanes + NumSegs + Offset + OffsetDist + 
                            OneWay + PaveType + PedCross + RTChannel + 
                            RTLanes + RTLnLength + RTMoveCtrl + RTWidth + 
                            Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                            TotalAADT + TotalT5YearInM + TurnProhib")

fmla.Offset <- as.formula("X5YrCrashCount ~ AADT + AreaType + City + County + IntAADT +
                            IntCat + IntTCType + LegRtType + LegSpeed + LegTCType + LegType + 
                            LegWidth + Lighting + LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                            MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                            NumLanes + NumSegs + Offset + OffsetDist + 
                            OneWay + PaveType + PedCross + RTChannel + 
                            RTLanes + RTLnLength + RTMoveCtrl + RTWidth + 
                            Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                            TotalT5YearInM + TurnProhib + offset(log(TotalAADT))") # offset(log())


#-- Poisson Regression
