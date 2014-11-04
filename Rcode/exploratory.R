################################################################################
# Exploratory Analysis
# To Fix: 
# - There are some predictors that are irrelevant, e.g.
  # in glm: problems since NumberLegs is given by IntCat
  # table(a$NumberLegs,a$IntCat)
# - In Exploratory Analysis:
#   - find importance of individual predictors (done in component.plots())
#   - find predictors with identical info
#     -- e.g. how well can x1 be predicted from x2, and x2 predicted from x1
#     -- this is something like mutual info
################################################################################

#------------------------------------------------
# Load Data and Packages
#------------------------------------------------
#-- Required Packages and Functions
library(dplyr)
library(spatstat)
library(ggmap) 
library(maptools)
library(scales)
library(mgcv)
source("Rcode/functions.R")  # has component.plots() function


#-- Load Accident Data 
load("data/accidents.RData") # loads a

#------------------------------------------------
# Make Intersection Data
# - data.frame x:
#   accidents - 5 year accident count at intersection
#   traffic - average traffic in 5 years (in millions)
#   rate - the accident rate per million vehicles
#   baseline - the expected number of accidents using overall state average
#   long/lat - spatial coordinates
# - avg.rate is average accident rate (over all collected data)
# Note: added new variable Traffic, which should replace TotalT5YeearInM
#       redo analysis below with Traffic replacing
#------------------------------------------------
library(dplyr)

#avg.rate = sum(a$X5YrCrashCount)/sum(a$TotalT5YearInM)

avg.rate = sum(a$X5YrCrashCount)/sum(a$Traffic) 

x = a %>% group_by(IntID) %>% summarize(
  accidents=sum(X5YrCrashCount),
  traffic=sum(Traffic), #sum(TotalT5YearInM),
  rate=accidents/traffic,
  baseline=avg.rate*traffic,
  Long=unique(Long),Lat=unique(Lat))


#------------------------------------------------
# Some basic plots of accidents vs. traffic vs. rate
# - note: lots of zeros in accidents, so be careful with taking logs
#------------------------------------------------
plot(x$traffic,x$accidents); abline(h=mean(a$X5YrCrashCount))
plot(x$traffic,log(x$rate)); abline(h=log(avg.rate))
plot(x$traffic,x$rate); abline(h=avg.rate)
plot(log(x$rate/avg.rate)); abline(h=0)
plot(x$baseline,log(.1+x$rate/avg.rate)); abline(h=0)
plot(x$baseline,x$accidents-x$baseline); abline(h=0)
plot(x$baseline,(x$accidents-x$baseline)/x$baseline); abline(h=0)

plot(x$traffic,log(.1+x$rate/avg.rate)); abline(h=0)
plot(x$traffic,x$accidents-x$baseline); abline(h=0)

"Total Traffic (in M)"
"Actual - Expected Accidents"
"Intersection / Overall Rate"



p1 = qplot(traffic,accidents,data=x) + 
    xlab("Total Traffic (in M)") + ylab("Number of Accidents") + geom_hline(yintercept=avg.rate)

p2 = qplot(traffic,log(.1+x$rate/avg.rate),data=x) + xlab("Total Traffic (in M)") +
    scale_size_area() +ylab("(Intersection / Overall) Rate") + geom_hline(yintercept=0)

p3 = qplot(traffic,accidents-baseline,data=x) + xlab("Total Traffic (in M)") +
    ylab("(Actual - Expected) Accidents") + geom_hline(yintercept=0)


#-- Our calculation from AADT doesn't relate well with accidents. IntAADT does
#   much better. Probably should use this for offset.

plot(a$TotalT5YearInM,a$X5YrCrashCount)
abline(lm(a$X5YrCrashCount~a$TotalT5YearInM),col=2)

plot(a$Traffic,a$X5YrCrashCount)
abline(lm(a$X5YrCrashCount~a$Traffic),col=2)




#-- Traffic rate vs. Accident Rate
#  Notice: clear distinction between major and minor roads, and accidents increase
#          as traffic increases (but for each type)
major = (a$LegType == 2)
minor = (a$LegType == 1)
plot(log(a$TotalT5YearInM),a$X5YrCrashCount,col=ifelse(major,"red","blue"))
points(log(a$TotalT5YearInM),loess(a$X5YrCrashCount~log(a$TotalT5YearInM),)$fitted,pch=19,cex=.5)
scatter.smooth(log(a$TotalT5YearInM),a$X5YrCrashCount,col=ifelse(major,"red","blue"))

scatter.smooth(log(a$TotalT5YearInM[major]),a$X5YrCrashCount[major],col="red")

scatter.smooth(log(a$TotalT5YearInM[minor]),a$X5YrCrashCount[minor],col="red")



#-- ggmap
library(ggmap)
alabama <- get_map("alabama", zoom = 7)
AlabamaMap <- ggmap(alabama, extent = "device", legend = 'topleft', inherit.aes = FALSE)

AlabamaMap + 
  geom_point(data=x, aes(x=Long, y=Lat, fill = log(accidents), size=log(rate)),
             shape=21,alpha=.8) +
  scale_fill_gradient(low="deepskyblue", high="red")
#  scale_fill_gradient2(low="deepskyblue",mid='black', high="red"),
#                       #midpoint=mean(x$accidents))


AlabamaMap + 
  geom_point(data=x, aes(x=Long, y=Lat, fill = accidents-baseline, size=log(traffic)),
             shape=21,alpha=.8) +
  scale_fill_gradient2(low="deepskyblue",mid='white', high="red",midpoint=0)


#------------------------------------------------
# Spatially smoothed plots of accident rate, accidents, and traffic
# - todo: 
#   - would be better to convert to equal area projection for smoothing
#   - get Rate Smooth Colors to have midpoint at avg.rate-or- convert rate plot 
#     to percent differnce from average
#------------------------------------------------
#-- use spatstat tools
library(spatstat)
library(maptools)
library(scales)
library(dplyr)

#-- Get Alabama Polygon
AL.map = map('state',"Alabama",fill=TRUE,plot=FALSE)
AL.poly = map2SpatialPolygons(AL.map, IDs=AL.map$names)
#library(rgdal)
#proj4string(AL.poly) = CRS("+proj=longlat")
#AL.poly = spTransform(AL.poly, CRS("+proj=utm +zone=16 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
AL = as.owin.SpatialPolygons(AL.poly)

#-- Make point pattern objects and Smoothing
pp = ppp(x$Long,x$Lat,window=AL,marks=dplyr::select(x,accidents,traffic,rate))
pp.accidents = ppp(x$Long,x$Lat,window=AL,marks=x$accidents)
pp.traffic = ppp(x$Long,x$Lat,window=AL,marks=x$traffic)

bw = .25
sm.accidents = Smooth(pp.accidents,sigma=bw,edge=FALSE)
sm.traffic = Smooth(pp.traffic,sigma=bw,edge=FALSE)
sm.rate = eval.im(sm.accidents/sm.traffic)


#-- Spatial Plots
get.colors = function(x) seq_gradient_pal(low="white",high="red")(rescale(x))
# see gradient_n_pal
# my.colors = colorRampPalette( c("blue","gray", "red"))

plot(sm.rate,main="Accident Rate (per million vehicles)")
plot(AL,add=TRUE)
points(pp,pch=19,col=get.colors(log(x$rate+.5)))

plot(sm.accidents,main='Accident Counts')
points(pp.accidents,pch=19,col=get.colors(log(marks(pp.accidents)+.5)),cex=.8)

plot(sm.traffic,main="Traffic Counts")
points(pp.traffic,pch=19,col=get.colors(marks(pp.traffic)))


#------------------------------------------------
# Sandbox: Get spatstat im to work with ggmap
#------------------------------------------------
# ?is there a way to add image to ggmap?
# https://groups.google.com/forum/#!topic/ggplot2/Bb8mq2emuSU
# http://www.obscureanalytics.com/2012/12/07/visualizing-baltimore-with-r-and-ggplot2-crime-data/
# https://stat.ethz.ch/pipermail/r-sig-geo/2014-April/020943.html

#-- Google Map with AL polygon shaded 
#   Notice that borders not exact! Google uses some sort of projection?
data = fortify(AL.poly)
AlabamaMap + geom_polygon(aes(x = long, y = lat), data = data,
                          color = "red", fill = "black", alpha = .4,size=1)  + theme_nothing()



#------------------------------------------------
# Component Plots - Find important variables 
# - todo: 
#   - try fam = negbin(c(1,10))
#   - test with offsets
#   - Need to get fmla from model_selection.R
#------------------------------------------------
library(mgcv)
source("Rcode/functions.R")  # has component.plots() function


# fmla <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
#                     LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
#                     LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
#                     MedWidth + MergeLanes + NextPIDist + NumberLegs + 
#                     NumLanes + Offset + OffsetDist + OneWay + PaveType + 
#                     PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
#                     RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
#                     log(TotalAADT) + log(TotalT5YearInM) + log(Traffic) +
#                     TurnProhib + Lat + Long") 
# 
# fmla.offset <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
#                     LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
#                     LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
#                     MedWidth + MergeLanes + NextPIDist + NumberLegs + 
#                     NumLanes + Offset + OffsetDist + OneWay + PaveType + 
#                     PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
#                     RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
#                     log(TotalAADT) + log(TotalT5YearInM) + offset(log(Traffic)) +
#                     TurnProhib + Lat + Long") 

fmla <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
                    LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
                    LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                    MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                    NumLanes + Offset + OffsetDist + OneWay + PaveType + 
                    PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
                    RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                    log(Traffic) + TurnProhib + Lat + Long") 

fmla.offset <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
                    LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
                    LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                    MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                    NumLanes + Offset + OffsetDist + OneWay + PaveType + 
                    PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
                    RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                    offset(log(Traffic)) + TurnProhib + Lat + Long") 

#- Variable Importance and Component Plots
# powerpoint = 1024 x 768
png("graph/componentplot-fmla.png",width=1024,height=700,units="px")
  score = component.plots(fmla,data=a)
dev.off()

png("graph/componentplot-offset.png",width=1024,height=700,units="px")
  score.offset = component.plots(fmla.offset,data=a)
dev.off()


varImportance = arrange(score,desc(score))
varI.offset = arrange(score.offset,desc(score))


png("graph/varImportance-barplot.png",width=1024,height=380,units="px")
#  par(mfrow=c(1,1),mar=c(4.8,7,3.5,.5))
#   barplot(score$score,names.arg=score$variable,
#           horiz=TRUE,las=1,cex.names=.75,
#           ylab="",xlab="Variable Importance",
#           main="No Offset: Componentwise Variable Importance (based on AIC)")
  par(mfrow=c(1,1),mar=c(7,4,2.5,.5))
  barplot(score$score,names.arg=score$variable,
          las=2,cex.names=.85,
          xlab="",ylab="Variable Importance",
          main="No Offset: Componentwise Variable Importance (based on AIC)")
dev.off()

png("graph/varImportance-barplot-offset.png",width=1024,height=380,units="px")
#   par(mfrow=c(1,1),mar=c(4.8,7,3.5,.5))
#   barplot(score.offset$score,names.arg=score.offset$variable,
#           horiz=TRUE,las=1,cex.names=.85,
#           ylab="",xlab="Variable Importance",
#           main="Offset: Componentwise Variable Importance (based on AIC)")
  par(mfrow=c(1,1),mar=c(7,4,2.5,.5))
  barplot(score.offset$score,names.arg=score.offset$variable,
          las=2,cex.names=.85,
          xlab="",ylab="Variable Importance",
          main="Offset: Componentwise Variable Importance (based on AIC)")
dev.off()




#----------------------------------------------------------