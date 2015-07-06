##########################################################
#  This script is to combine and clean raw csv files into
#  a .rdata file.
#
#  Author: Ketong Wang
##########################################################

#-- Read in the data
a <- read.csv("data/accident.csv")
vClass <- read.csv("data/classInfo.csv")
vClass[vClass$VarName=="MergeLanes","Class"] = "categorical"
coordinate <- read.csv("data/coordinate.csv")
a[(a == "<Null>")] = NA  # convert <Null> to NA

#-- Remove data from IntCat==9
a = subset(a,!(IntCat == 9))  # This should remove all the intersections and legs

#-- convert categorical data to "factor"
ind = which(vClass$Class=="categorical")
for(j in ind) a[,j] = factor(a[,j])

#-- convert numeric data to "numeric
ind = which(vClass$Class=="numeric")
for(j in ind) class(a[,j]) = "numeric"

#-- replace missing AADT
missing.AADT = is.na(a$AADT)   # indicates if AADT is missing (NA)
a$AADT[missing.AADT] = a$IntAADT[missing.AADT]

#-- Correct NumberLegs data
a$NumberLegs[a$NumberLegs == 44] = 4

#-- Create intersection related variables other intersection 
#   variables to consider: variation in AADT at intersection, 
#   total volume at intersection, etc

b = strsplit(as.character(a$LegID),"_")
b = do.call(rbind,b)
IntID = factor(b[,1])    # intersection ID variable
segTab = as.data.frame(table(IntID))
ind = match(IntID,segTab$IntID)
NumSegs = as.numeric(segTab$Freq[ind])
total = tapply(a$AADT,IntID,sum)  
TotalAADT = as.numeric(total[ind]) # total traffic at intersection

a = data.frame(IntID,a,TotalAADT)  # add new variables

#-- One more variable leg 5 year total traffic per million
#   Note: AADT is not related well to accidents, better measure is IntAADT
a$TotalT5YearInM <- a$AADT * 5 * 365 / 1000000
a$Traffic <- a$IntAADT * 5 * 365 / 1000000  # this is much better related to accidents

#-- Merge and save accident data with coordinate data
a <- merge(a, coordinate, by = "IntID")

write.csv(a,"data/cleaned.csv", row.names=FALSE)
save(a, file = "data/accidents.RData")


