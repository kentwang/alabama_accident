################################################################################
# Exploratory Analysis
# - 
################################################################################
#setwd("alabama_accident")


#------------------------------------------------
# Load Data and Packages
#------------------------------------------------
#-- Required Packages and Functions
library(dplyr)
library(spatstat)
library(ggmap)
library(maptools)
library(scales)

#-- Load Accident Data 
load("data/accidents.RData") # loads a

#  a$TotalT5YearInM <- a$AADT * 5 * 365 / 1000000


#------------------------------------------------
# Make Intersection Data
# - data.frame x:
#   accidents - 5 year accident count at intersection
#   traffic - average traffic in 5 years (in millions)
#   rate - the accident rate per million vehicles
#   baseline - the expected number of accidents using overall state average
#   long/lat - spatial coordinates
# - avg.rate is average accident rate (over all collected data)
#------------------------------------------------
library(dplyr)

avg.rate = sum(a$X5YrCrashCount)/sum(a$TotalT5YearInM)

x = a %>% group_by(IntID) %>% summarize(
          accidents=sum(X5YrCrashCount),
          traffic=sum(TotalT5YearInM),
          rate=accidents/traffic,
          baseline=avg.rate*traffic,
          Long=unique(Long),Lat=unique(Lat))


#------------------------------------------------
# Some basic plots of accidents vs. traffic vs. rate
# - note: lots of zeros in accidents, so be careful with taking logs
#------------------------------------------------
plot(x$traffic,x$accidents); abline(h=mean(a$X5YrCrashCount))
plot(x$traffic,log(x$rate+.01)); abline(h=log(avg.rate))
plot(x$traffic,x$rate); abline(h=avg.rate)
plot(log(x$rate/avg.rate)); abline(h=0)
plot(x$baseline,log(x$rate/avg.rate)); abline(h=0)
plot(x$baseline,x$accidents-x$baseline); abline(h=0)
plot(x$baseline,(x$accidents-x$baseline)/x$baseline); abline(h=0)

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

#-- Get Alabama Polygon
AL.map = map('state',"Alabama",fill=TRUE,plot=FALSE)
AL.poly = map2SpatialPolygons(AL.map, IDs=AL.map$names)
  #library(rgdal)
  #proj4string(AL.poly) = CRS("+proj=longlat")
  #AL.poly = spTransform(AL.poly, CRS("+proj=utm +zone=16 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
AL = as.owin.SpatialPolygons(AL.poly)

#-- Make point pattern objects and Smoothing
pp = ppp(x$Long,x$Lat,window=AL,marks=select(x,accidents,traffic,rate))
pp.accidents = ppp(x$Long,x$Lat,window=AL,marks=x$accidents)
pp.traffic = ppp(x$Long,x$Lat,window=AL,marks=x$traffic)

bw = .15
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








