#############################################
## A makefile for all the graphs in AAP paper
#############################################

library(dplyr)
library(spatstat)
library(ggmap) 
library(maptools)
library(scales)
library(mgcv)
library(glmnet)
library(gplots)
source("Rcode/functions.R")  # has component.plots() function
source('Rcode/glm_my.R')



#-- Load Accident Data 
load("data/accidents.RData") # loads a


#======================================#
#-- There is a potential outlier in observation 63
  plot(log(a$Traffic),a$X5YrCrashCount)
  # one observation has 49 accidents, all others <= 21! 
  # is this an outlier, or just bad intersection?
  plot(log(a$Traffic),log(a$X5YrCrashCount)) # doesn't *look* so bad, but outlier will be big in likelihood

#- Should it be removed, or is there something special about this intersection leg?
# outlier = which(a$X5YrCrashCount>21)
# a = a[-outlier,]

data = model.frame(fmla.offset,data=a)  
offset = as.vector(model.offset(data))
y = model.response(data)    # response vector
vars = attr(terms(data),"term.labels")
par(mfrow=c(7,6))
for(j in 1:length(vars)){
  var = vars[j]
  plot(data[,var],y,xlab=paste(var))
  points(data[outlier,var],y[outlier],col=2,cex=2,pch=19)
}
    
#======================================#


#-- Create formula
fmla.offset <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
                    LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
                    LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                    MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                    NumLanes + Offset + OffsetDist + OneWay + PaveType + 
                    PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
                    RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                    offset(log(Traffic)) + TurnProhib + Lat + Long") 


#--------------------------------------------
# Figure 1 Exploratory Analysis
#--------------------------------------------


  data = model.frame(fmla.offset,a)  # make correct transformations from fmla
  offset = as.vector(model.offset(data))
  y = model.response(data)    # response vector
  vars = attr(terms(data),"term.labels") # predictor variables
  p = length(vars)


#-- find component scores
score = edf = numeric(p)
for(j in 1:p){
  x = data[,vars[j]]
  fit = component.fit(x,y,offset,fam=poisson(),max.df=6,plot=FALSE)
  score[j] = fit$score
  edf[j] = fit$edf
}
score = data.frame(vars,score,edf,var.type=sapply(data[,vars],class),
                  stringsAsFactors=FALSE)  
ordered_component = score$vars[order(-score$score)]  # decreasing order of score

score[ordered_component,]   
  
  
#-- plot all component plots (ordered from largest to lowest score)
par(mfrow=c(7,6),mar=c(2,2,1.5,1))
for(j in 1:p){
  var = ordered_component[j]
  x = data[,var]
  fit = component.fit(x,y,offset,fam=poisson(),max.df=6,plot=TRUE)
  title(paste0(var," (",round(fit$score),")"))
}




#-- Draw four most and two least influential factors
VARS = c("IntCat", "LegRtType", "SkewAngle", "RTLnLength")
#score.offset = component.plots(fmla.offset,data=a,VARS=VARS,xy=c(2,2))


data = model.frame(fmla.offset,data=a)  # make correct transformations from fmla
offset = as.vector(model.offset(data))
y = model.response(data)    # response vector
vars = attr(terms(data),"term.labels")

par(mfrow=c(2,2),mar=c(4,4,1.5,1))
for(i in 1:length(VARS)){
  var = VARS[i]
  fit = component.fit(data[,var],y,offset=offset,fam=poisson(),max.df=6, #min.sp=500,
                            plot=TRUE, ylim=c(-2,2)) 
  mtext(side = 1, paste0("(",letters[i], ") ",var," (",round(fit$score),")"), 
        line = 2.5, cex = 1.2)
}



           
              
#--------------------------------------------
# Figure 2 -mlogL predictive evaluation
#--------------------------------------------
fold = cvfolds(nrow(a),k=20,seed=9122014)  # get cv partition info fold 20 and 5
Y = a$X5YrCrashCount  # response

max.complex = 10
tree.seq = seq(500,10000,by=100)
cp.seq = seq(0, 0.05, by = 0.0001)

# fit models with cross-validation
mu.cv.poisReg.offset <- cv.poisReg(fmla.offset, data = a, fold= fold, max.complex = max.complex)
mu.cv.negBino.offset <- cv.negBino(fmla.offset, data = a, fold= fold, max.complex = max.complex)
mu.cv.glmnet.offset = cv.glmnetPois(fmla.offset,data=a,fold=fold)
mu.cv.tree.offset = cv.treePois(fmla.offset, data = a, fold = fold, cp.seq = cp.seq, show.pb=TRUE)
set.seed(20150702)
gbm_4.offset <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                           tree.seq = tree.seq, interaction.depth=4,
                           shrinkage = .005)

# complexity tuning parameters
cp.poisReg = 1:max.complex
cp.negBino = 1:max.complex
cp.glmnet = attr(mu.cv.glmnet.offset, "lambda")
cp.glmnet = -log(cp.glmnet/max(cp.glmnet))
#cp.glmnet = cp.glmnet + (10-max(cp.glmnet))
cp.glmnet = 10*cp.glmnet/max(cp.glmnet)
cp.gbm = tree.seq/1000
cp.tree = cp.seq*200

# Figure 2

model.info = data.frame(
  model= c("Poisson", "NB", "glmnet", "Trees", "BRT"),
  pch = c(15,16,NA_integer_,NA_integer_,NA_integer_),
  col = c(6, 2, 3, 'orange', 4),
  lty = c(1,1,1,4,2),
  pt.cex = c(1.2,1.2,1,1,1),stringsAsFactors=FALSE
  )



par(mfrow=c(1,1),mar=c(4.1,4,1,1))
plot(c(1,max.complex),c(1.5,2.1),typ='n',
     ylab="prediction error",xlab="model complexity")
with(model.info[1,], lines(1:max.complex, mlogL(mu.cv.poisReg.offset, Y),
                           type='o',pch=pch,cex=pt.cex,col=col,lty=lty))
with(model.info[2,], lines(1:max.complex, mlogL(mu.cv.negBino.offset, Y),
                           type='o',pch=pch,cex=pt.cex,col=col,lty=lty))
with(model.info[3,], lines(cp.glmnet, na.omit(mlogL(mu.cv.glmnet.offset, Y)),
                           type='o',pch=pch,cex=pt.cex,col=col,lty=lty))
with(model.info[5,], lines(cp.gbm, na.omit(mlogL(gbm_4.offset, Y)),
                           type='o',pch=pch,cex=pt.cex,col=col,lty=lty))
with(model.info[c(1,2,3,5),],
legend("topleft" ,legend=model, col=col,pch=pch,lty=lty,pt.cex=pt.cex,
       lwd = 2, ncol = 1, cex = 1.2,bty='n') #text.font = 3
)

#plot(1:max.complex, mlogL(mu.cv.poisReg.offset, Y), type = 'b', pch = 15, cex = 1.2, 
#     col=6,lwd=2, ylab="prediction error", ylim = c(1.57, 2.2),  xlab="model complexity",
#     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.2)
#lines(1:max.complex, mlogL(mu.cv.negBino.offset, Y), type = 'b', pch = 16, cex = 1.2, col=2,lwd=2)
#lines(cp.glmnet, na.omit(mlogL(mu.cv.glmnet.offset, Y)),col=3,lwd=2)
#lines(cp.gbm,mlogL(gbm_4.offset,Y),col=4,lty=2,lwd=2)
# lines(cp.tree, mlogL(mu.cv.tree.offset, Y), col='orange',lty=4,lwd=2)
# legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "Trees", "BRT"), 
#        col = c(6, 2, 3, 'orange', 4), pch=c(15,16,NA_integer_,NA_integer_,NA_integer_),
#        lty=c(1,1,1,4,2),lwd = 2, text.font = 3, ncol = 5, pt.cex=c(1.2,1.2,1,1,1), cex = 1.2)




#--------------------------------------------
# Table 2. Optimal complexities
#--------------------------------------------
cp.poisReg = (1:max.complex)[which.min( mlogL(mu.cv.poisReg.offset , Y))]
cp.negBino = (1:max.complex)[which.min( mlogL(mu.cv.negBino.offset , Y))]
cp.glmnet = (attr(mu.cv.glmnet.offset, "lambda"))[which.min(mlogL(mu.cv.glmnet.offset , Y))]
cp.gbm = tree.seq[which.min(mlogL(gbm_4.offset , Y))]
cp.tree = cp.seq[which.min(mlogL(mu.cv.tree.offset , Y))]


#--------------------------------------------
# Figure 3. Factor importance of the best models
#--------------------------------------------
# fit best models
best.poisReg = best.poisReg.fit(fmla.offset, data = a, cp.poisReg = cp.poisReg)
best.negBino = best.negBino.fit(fmla.offset, data = a, cp.negBino = cp.negBino)

X = model.matrix(fmla.offset,a)[,-1]
Y = as.matrix(a[,as.character(fmla.offset[[2]])])
offset = model.offset(model.frame(fmla.offset,a))
best.glmnet = glmnet(X,Y,offset=offset,family="poisson",alpha=0.8, lambda = cp.glmnet)

best.gbm = summary(gbm(fmla.offset, data=a, distribution = "poisson", n.trees = cp.gbm,
                       interaction.depth = 4), plotit=FALSE)

tree.0 = rpart(fmla.offset, data = a, method = "poisson", cp = 0,xval=0, minbucket=3)
best.tree = prune(tree.0,cp=cp.tree)

###### Cleveland dot plot
fmla.string = strsplit(as.character(fmla.offset)[3], " \\+ ")[[1]]
fmla.string = fmla.string[fmla.string != "offset(log(Traffic))"]
allVariables = fmla.string[order(fmla.string)]

# reorder the importance using the order of the component graph
ordered_component = c("IntCat", "LegRtType", "AreaType", "LegWidth", "OneWay", "Lat",
                      "Terrain", "Long", "NextPIDist", "Lighting", "LegTCType", "LegSpeed",
                      "MedWidth", "LegType", "OffsetDist", "RTMoveCtrl", "Offset",
                      "LTWidth", "IntTCType", "LTOffset", "NumLanes", "SightRt",
                      "LTLanes", "LTLnLength", "PedCross", "RTWidth", "SightLt",
                      "PaveType", "TurnProhib", "MedType", "RTLanes", "RTChannel",
                      "MergeLanes", "NumberLegs", "Rumble", "SkewAngle", "RTLnLength")

comp_order = match(ordered_component, allVariables) # indices from allvariabels to ordered_component

varImp <- as.vector(c(varImpStandard(summary(best.poisReg)$coefficients[-1, 3], allVariables)[comp_order],
                      varImpStandard(summary(best.negBino)$coefficients[-1, 3], allVariables)[comp_order],
                      varImpStandard(best.glmnet.imp ,allVariables)[comp_order],
                      varImpStandard(best.tree$variable.importance, allVariables)[comp_order],
                      varImpStandard(best.gbm.inf, allVariables)[comp_order]))

model.names <- c("Poisson", "NB",
                 "glmnet", "Trees", "BRT")

importances <- as.data.frame(cbind(rep(ordered_component, length(model.names)),
                                   rep(model.names, each = length(allVariables))))

importances$importance <- varImp
names(importances) <- c("var", "model", "importance")

# Cleveland Dot Plot using ggplot
dot_theme = theme_bw() +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(color="grey60",
                                        linetype="dashed"))


importances_trunc = importances[!importances$var %in% c("MedWidth", "Lighting", "LegSpeed", "SkewAngle", "LegType", "NumLanes", "OffsetDist", "SightLt", "SightRt", "RTLnLength", "RTWidth", "LTOffset", "RTChannel","RTLanes", "Rumble", "PaveType", "MergeLanes"), ]

ggplot(importances_trunc, aes(y=reorder(var, importance, max), x=importance)) + 
  geom_point(aes(shape = model, col = model), size = 3) + 
  dot_theme +
  theme(legend.position = c(0.85, .77), legend.background = element_rect(colour = "black")) +
  scale_y_discrete("")


#--------------------------------------------
# Figure 4. Interaction effect
#--------------------------------------------
optimal.gbm = gbm(fmla.offset, data = a, distribution = "poisson", n.trees = 3700, interaction.depth = 4)
vars = attr(terms(model.frame(fmla.offset, a)),"term.labels")
p = length(vars)

interaction.mag = matrix(0, p, p)

for(i in 1:(p-1)) {
  for(j in (i+1):p) {
    cat(i, j, "\n")
    interaction.mag[i, j] = interact.gbm(optimal.gbm, data = a, i.var = c(i, j), n.trees = 3700)
  }
}

colnames(interaction.mag) = vars
rownames(interaction.mag) = vars

interaction.mag[is.na(interaction.mag)] = 0


BRT.vars = c("Terrain", "SkewAngle", "PedCross", "NextPIDist", "MedWidth", 
             "Long", "LegWidth", "LegRtType", "Lat", "IntCat") # from the BRT variable importances

BRT.vars = c("Terrain", "SkewAngle", "PedCross", "NextPIDist", "MedWidth", "LegWidth", "LegRtType", "IntCat")
interaction.index = which(vars %in% BRT.vars)
BRT.vars = vars[interaction.index]

my_palette <- colorRampPalette(c("white", "blue"))(n = 299)

par(cex.main=0.5)

# save as portrait 12 by 6
heatmap.2(interaction.mag[BRT.vars, BRT.vars], dendrogram = "none", Colv = FALSE,
          Rowv = FALSE, key = T, keysize = 1, density.info="none", trace="none",
          key.xlab = "Scale of H-Statistic", key.title = NA,
          cexRow = 1.4, cexCol = 1.4, margins = c(7, 9), col = my_palette,
          srtCol=30, main = "Two-way factor interactions")











