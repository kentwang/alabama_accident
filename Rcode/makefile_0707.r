#############################################
## A makefile for all the graphs in AAP paper
#############################################

library(ggplot2)
library(tidyr)
library(dplyr)
#library(scales)
library(mgcv)
library(glmnet)
library(gplots)
source("Rcode/functions.R")  # has component.plots() function
source('Rcode/glm_my.R')




#-- Plotting Properties
plotDir = 'graph_nolatlong'
# plotDir = 'graph'
#plotDir = 'C:/GoogleDrive/PROJECTS/Accident/paper/paper_latex/graph'

pdf.options(family='Bookman',pointsize=10)
setEPS(horizontal=FALSE,paper='special',family='Bookman',
       pointsize=12,width=10,height=6)



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


a %>% group_by(IntID) %>% summarize(avg=mean(X5YrCrashCount)) %>% arrange(desc(avg))
a %>% filter(IntID %in% c('414350','125862'))
#outlier = which(a$IntID %in% c('414350','125862'))
#a = a[-outlier,]

#======================================#


#-- Create formula
# fmla.offset <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
#                     LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
#                     LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
#                     MedWidth + MergeLanes + NextPIDist + NumberLegs + 
#                     NumLanes + Offset + OffsetDist + OneWay + PaveType + 
#                     PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
#                     RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
#                     offset(log(Traffic)) + TurnProhib + Lat + Long") 

fmla.offset <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
                    LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
                    LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                    MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                    NumLanes + Offset + OffsetDist + OneWay + PaveType + 
                    PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
                    RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                    offset(log(Traffic)) + TurnProhib") # formula without lat/long



#--------------------------------------------
# Exploratory Analysis
#--------------------------------------------

  data = model.frame(fmla.offset,a)  # make correct transformations from fmla
  offset = as.vector(model.offset(data))
  y = model.response(data)    # response vector
  vars = attr(terms(data),"term.labels") # predictor variables
  p = length(vars)


#-- find component scores
score = edf = nlevs = numeric(p)
for(j in 1:p){
  x = data[,vars[j]]
  fit = component.fit(x,y,offset,fam=poisson(),max.df=6,plot=FALSE)
  score[j] = fit$score           
  edf[j] = fit$edf               # effective degrees of freedom 
#  nlevs[j] = nlevels(x)          
  nlevs[j] = ifelse(is.factor(x), # number of unique (or possible) values
                    nlevels(x),
                    length(unique(x)))
  
}
score = data.frame(vars,score,edf,var.type=sapply(data[,vars],class),nlevs,
                  stringsAsFactors=FALSE)  
ordered_component = score$vars[order(-score$score)]  # decreasing order of score




#--------------------------------------------
# Basic Info
#--------------------------------------------
nrow(a)                  # number of legs
length(unique(a$IntID))  # number of intersections
length(attr(terms(fmla.offset),"term.labels"))  # number of predictors
length(attr(terms(fmla.offset),"offset"))       # offset indicator


#--------------------------------------------
# Table 1 (old) and new separate tables
#--------------------------------------------
transform(score[ordered_component,],score=round(score),edf=round(edf,1)) # old Table 1  

original_table = read.table("data/original_table_1_nobackslash.txt", sep='&')
row.names(original_table) = original_table$V1

cat('Factor Name & Definition & Level & Strength & Min & Max & Mean & SD', "\n")
for(i in 1:length(ordered_component)) {
  var = ordered_component[i]
  varCol = data[, var]
  if(!is.factor(varCol)) {
    cat('\\textit{', trimws(as.character(original_table[var, 'V1'])), '}', '&', 
    trimws(as.character(original_table[var, 'V2'])), '&',
    trimws(as.character(original_table[var, 'V3'])), '&',
    score[var, 'score'],  "&", min(varCol), "&", max(varCol), "&", 
    round(mean(varCol), 2), "&", round(sd(varCol), 2), "\\\\\n", sep='')
  }
}

cat('Factor Name & Definition & Level & Strength & \\# categories', "\n")
for(i in 1:length(ordered_component)) {
  var = ordered_component[i]
  varCol = data[, var]
  if(is.factor(varCol)) {
    cat('\\textit{', trimws(as.character(original_table[var, 'V1'])), '}', '&', 
    trimws(as.character(original_table[var, 'V2'])), '&',
    trimws(as.character(original_table[var, 'V3'])), '&', 
    score[var, 'score'], "&", nlevels(varCol), "\\\\\n", sep='')
  }
}



#--------------------------------------------
# Figure 1 Exploratory Analysis
#--------------------------------------------
#-- Draw four most and two least influential factors
VARS = c("IntCat", "OneWay", "SkewAngle", "RTLnLength")
#score.offset = component.plots(fmla.offset,data=a,VARS=VARS,xy=c(2,2))


data = model.frame(fmla.offset,data=a)  # make correct transformations from fmla
offset = as.vector(model.offset(data))
y = model.response(data)    # response vector
vars = attr(terms(data),"term.labels")



pdf(file=file.path(plotDir,"Fig1.pdf"),width=8,height=6)
# jpeg(file=file.path(plotDir,"Fig1.jpeg"),width=800,height=600,quality=150)
  par(mfrow=c(2,2),mar=c(4,2,1.5,1))
  for(i in 1:length(VARS)){
    var = VARS[i]
    fit = component.fit(data[,var],y,offset=offset,fam=poisson(),max.df=6, #min.sp=500,
                              plot=TRUE, ylim=c(-2,2),ylab="") 
    mtext(side = 1, paste0("(",letters[i], ") ",var," (",round(fit$score),")"), 
          line = 2.5, cex = 1.2)
  }
dev.off()




#-- plot all component plots (ordered from largest to lowest score)
pdf(file=file.path(plotDir,"Complete_Component.pdf"),width=8,height=6)
par(mfrow=c(7,6),mar=c(2,2,1.5,1))
for(j in 1:p){
  var = ordered_component[j]
  x = data[,var]
  fit = component.fit(x,y,offset,fam=poisson(),max.df=6,plot=TRUE)
  title(paste0(var," (",round(fit$score),")"))
}
dev.off()

              
#--------------------------------------------
# Figure 2 -mlogL predictive evaluation
#--------------------------------------------
fold = cvfolds(nrow(a),k=20,seed=9122014)  # get cv partition info fold 20 and 5
Y = a$X5YrCrashCount  # response

max.complex = 10
tree.seq = seq(300,3000,by=100)
cp.seq = seq(0, 0.05, by = 0.0001)

# fit models with cross-validation
mu.cv.poisReg.offset <- cv.poisReg(fmla.offset, data = a, fold= fold, max.complex = max.complex)
mu.cv.negBino.offset <- cv.negBino(fmla.offset, data = a, fold= fold, max.complex = max.complex)
mu.cv.glmnet.offset = cv.glmnetPois(fmla.offset,data=a,fold=fold)
#mu.cv.tree.offset = cv.treePois(fmla.offset, data = a, fold = fold, cp.seq = cp.seq, show.pb=TRUE)
set.seed(20150702)
idepth = 7           # ensure same settings for best.gbm (line 341)
shrink = .005
gbm_7.offset_01  <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                           tree.seq = tree.seq, interaction.depth=idepth,
                           shrinkage = shrink)
gbm_4.offset <- gbm_7.offset_01 # use this assignment without changing notation

# complexity tuning parameters
complex.poisReg = 1:max.complex
complex.negBino = 1:max.complex
complex.glmnet = -1.2*log(attr(mu.cv.glmnet.offset, "lambda"))
#complex.glmnet = -log(complex.glmnet/max(complex.glmnet))
#complex.glmnet = complex.glmnet + (10-max(complex.glmnet))
#complex.glmnet = 10*complex.glmnet/max(complex.glmnet)
complex.gbm = tree.seq/300
complex.tree = cp.seq*200

# Figure 2

plot.info = data.frame(
  model= c("Poisson", "NB", "glmnet", "BRT"),
  #pch = c(15,16,NA_integer_,NA_integer_,NA_integer_),
  pch = c(15,16,17,18),
  pt.cex = c(1.2,1.4,1.2,1.5),
  #col = c(6, 2, 3, 'orange', 4),
  #col = c("#668F3C","#BA57B0","#C5583F","#617DB4"),
  col = c("#e41a1c","#377eb8","#4daf4a","#984ea3"), # colorbrewer
  lty = c(1,1,1,1),
  stringsAsFactors=FALSE)
  



pdf(file=file.path(plotDir,"Fig2.pdf"),width=8,height=3)
  par(mfrow=c(1,1),mar=c(4.1,4,1,1),las=1)
  plot(c(1,max.complex),c(1.57,2.2),typ='n',
       ylab="prediction error",xlab="model complexity",xaxt='n')
  axis(1,1:10)
  with(plot.info[1,], lines(complex.poisReg, mlogL(mu.cv.poisReg.offset, Y),
                             type='o',pch=pch,cex=pt.cex,col=col,lty=lty))
  with(plot.info[2,], lines(complex.negBino, mlogL(mu.cv.negBino.offset, Y),
                             type='o',pch=pch,cex=pt.cex,col=col,lty=lty))
  
  pt.seq = unique(round(seq(1,length(complex.glmnet),length=11)))
  with(plot.info[3,], lines(complex.glmnet, na.omit(mlogL(mu.cv.glmnet.offset, Y)),
                             pch=pch,cex=pt.cex,col=col,lty=lty))
  with(plot.info[3,], points(complex.glmnet[pt.seq], na.omit(mlogL(mu.cv.glmnet.offset, Y))[pt.seq],
                             pch=pch,cex=pt.cex,col=col))
  
  pt.seq = unique(round(seq(1,length(complex.gbm),length=11)))
  with(plot.info[4,], lines(complex.gbm, na.omit(mlogL(gbm_4.offset, Y)),
                             pch=pch,cex=pt.cex,col=col,lty=lty))
  with(plot.info[4,], points(complex.gbm[pt.seq], na.omit(mlogL(gbm_4.offset, Y))[pt.seq],
                             pch=pch,cex=pt.cex,col=col))
  with(plot.info,
  legend("topleft" ,legend=model, col=col,pch=pch,lty=lty,pt.cex=pt.cex,
         lwd = 2, ncol = 1, cex=1,bty='n') #text.font = 3
  )
dev.off()



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

#-- mdp 8/31/2015 don't think these are needed
# score.poisReg = data.frame(complexity=1:max.complex,c2=1:max.complex,score=mlogL(mu.cv.poisReg.offset,Y))
# score.negBino = data.frame(complexity=1:max.complex,c2=1:max.complex,score=mlogL(mu.cv.negBino.offset,Y))
# score.glmnet = data_frame(complexity=attr(mu.cv.glmnet.offset, "lambda"),
#                           c2=-log(complexity/max(complexity)),
#                           score=na.omit(mlogL(mu.cv.glmnet.offset,Y)))
# ## Note: another step with complexity.
# # score.tree = data_frame(complexity=cp.seq,c2=cp.seq/200,score=mlogL(mu.cv.tree.offset, Y))
# score.gbm = data_frame(complexity=tree.seq,c2=tree.seq/1000,score=na.omit(mlogL(gbm_4.offset, Y)))
# 




#--------------------------------------------
# Table 2. Optimal complexities
#--------------------------------------------
cp.poisReg = (1:max.complex)[which.min( mlogL(mu.cv.poisReg.offset , Y))]
cp.negBino = (1:max.complex)[which.min( mlogL(mu.cv.negBino.offset , Y))]
cp.glmnet = (attr(mu.cv.glmnet.offset, "lambda"))[which.min(mlogL(mu.cv.glmnet.offset , Y))]
cp.gbm = tree.seq[which.min(mlogL(gbm_4.offset , Y))]
#cp.tree = cp.seq[which.min(mlogL(mu.cv.tree.offset , Y))]


data.frame(
 model = c("Poisson", "NB", "glmnet", "BRT"),
 transformed = c(cp.poisReg,cp.negBino,-1.2*log(cp.glmnet),cp.gbm/300),
 original = c(cp.poisReg,cp.negBino,cp.glmnet,cp.gbm),
 prediction.error = round( c(min( mlogL(mu.cv.poisReg.offset , Y)),
                      min( mlogL(mu.cv.negBino.offset , Y),na.rm=TRUE),
                      min(mlogL(mu.cv.glmnet.offset , Y),na.rm=TRUE),
                      min(mlogL(gbm_4.offset , Y),na.rm=TRUE) ),2)
)

   

#--------------------------------------------
# Figure 3. Factor importance of the best models using full data
#--------------------------------------------
# fit best models
best.poisReg = best.poisReg.fit(fmla.offset, data = a, cp.poisReg = cp.poisReg)
best.negBino = best.negBino.fit(fmla.offset, data = a, cp.negBino = cp.negBino)

X = model.matrix(fmla.offset,a)[,-1]
Y = as.matrix(a[,as.character(fmla.offset[[2]])])
offset = model.offset(model.frame(fmla.offset,a))
best.glmnet = glmnet(X,Y,offset=offset,family="poisson",alpha=0.8, lambda = cp.glmnet)

# best.gbm = summary(gbm(fmla.offset, data=a, distribution = "poisson", n.trees = cp.gbm,
#                        interaction.depth = 4), plotit=FALSE)
best.gbm = gbm(fmla.offset, data=a, distribution = "poisson", n.trees = cp.gbm,
                        interaction.depth = idepth, shrinkage=shrink)



#tree.0 = rpart(fmla.offset, data = a, method = "poisson", cp = 0,xval=0, minbucket=3)
#best.tree = prune(tree.0,cp=cp.tree)

#-- Added by Ketong --#
# processing of glmnet
# Todo standardization of coefficient by sd(X)
best.glmnet.imp = as.vector(best.glmnet$beta)/apply(X, 2, sd) # standardized by the coded model matrix
#names(best.glmnet.imp) = rownames(best.glmnet$beta)

# processing of gbm
best.gbm.inf = summary(best.gbm,plot=FALSE)$rel.inf
names(best.gbm.inf) = summary(best.gbm,plot=FALSE)$var
#-- Added by Ketong --#


################################

varimp <- function(x,model="glm",fmla=fmla.offset,data=a,varlist=NULL){
  if(model=="glm"){
    vars = attr(x$terms,"term.labels")
    z = abs(summary(x)$coefficients[,3])  # abs(z-value)
    v = data.frame(var=vars,
                   rel.inf=sapply(vars,function(x) max(z[grep(x,names(z))])))
  }
  if(model=="gbm"){
    v = summary(x,plot=FALSE)
  }
  if(model=="glmnet"){
    X = model.matrix(fmla,data)[,-1]
    vars = attr(terms(fmla),"term.labels")
    z = abs(x$beta)/apply(X,2,sd)
    v = data.frame(var=vars,
                   rel.inf=sapply(vars,function(x) max(z[grep(x,rownames(z))])))
  }
  v = transform(v,rel.inf=100*rel.inf/sum(rel.inf))  # sum to 100
  if(!is.null(varlist)){
    rel.inf = numeric(length(varlist))
    ind = match(v$var,varlist)
    rel.inf[ind] = v$rel.inf
    v = data.frame(var=varlist,rel.inf)
  }
return(v)
}


tmp2 = data.frame(
  var = ordered_component,
  Poisson = varimp(best.poisReg,model="glm",varlist=ordered_component)$rel.inf,
  NB = varimp(best.negBino,model="glm",varlist=ordered_component)$rel.inf,
  glmnet = varimp(best.glmnet,model="glmnet",fmla=fmla.offset,data=a,varlist=ordered_component)$rel.inf,
  BRT = varimp(best.gbm,model="gbm",varlist=ordered_component)$rel.inf
)
#######################################



###### Cleveland dot plot
#fmla.string = strsplit(as.character(fmla.offset)[3], " \\+ ")[[1]]
#fmla.string = fmla.string[fmla.string != "offset(log(Traffic))"]
#allVariables = fmla.string[order(fmla.string)]
#allVariables = attr(terms(fmla.offset),"term.labels")

#- This was set above in Table 1
# reorder the importance using the order of the component graph
# ordered_component = c("IntCat", "LegRtType", "AreaType", "LegWidth", "OneWay", "Lat",
#                       "Terrain", "Long", "NextPIDist", "Lighting", "LegTCType", "LegSpeed",
#                       "MedWidth", "LegType", "OffsetDist", "RTMoveCtrl", "Offset",
#                       "LTWidth", "IntTCType", "LTOffset", "NumLanes", "SightRt",
#                       "LTLanes", "LTLnLength", "PedCross", "RTWidth", "SightLt",
#                       "PaveType", "TurnProhib", "MedType", "RTLanes", "RTChannel",
#                       "MergeLanes", "NumberLegs", "Rumble", "SkewAngle", "RTLnLength")
#
# comp_order = match(ordered_component, allVariables) # indices from allvariabels to ordered_component

# varImp <- as.vector(c(varImpStandard(summary(best.poisReg)$coefficients[-1, 3], allVariables)[comp_order],
#                       varImpStandard(summary(best.negBino)$coefficients[-1, 3], allVariables)[comp_order],
#                       varImpStandard(best.glmnet.imp ,allVariables)[comp_order],
#                       varImpStandard(best.tree$variable.importance, allVariables)[comp_order],
#                       varImpStandard(best.gbm.inf, allVariables)[comp_order]))
# 
# model.names <- c("Poisson", "NB",
#                  "glmnet", "Trees", "BRT")
# 
# importances <- as.data.frame(cbind(rep(ordered_component, length(model.names)),
#                                    rep(model.names, each = length(allVariables))))
# 
# importances$importance <- varImp
# names(importances) <- c("var", "model", "importance")
# 
# # Cleveland Dot Plot using ggplot
dot_theme = theme_bw() +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_line(color="grey60",
                                        linetype="dashed"))
# 
# 
# importances_trunc = importances[!importances$var %in% c("MedWidth", "Lighting", "LegSpeed", "SkewAngle", "LegType", "NumLanes", "OffsetDist", "SightLt", "SightRt", "RTLnLength", "RTWidth", "LTOffset", "RTChannel","RTLanes", "Rumble", "PaveType", "MergeLanes"), ]
# 
# ggplot(importances_trunc, aes(y=reorder(var, importance, max), x=importance)) + 
#   geom_point(aes(shape = model, col = model), size = 3) + 
#   dot_theme +
#   theme(legend.position = c(0.85, .77), legend.background = element_rect(colour = "black")) +
#   scale_y_discrete("") 


##############################################################
varImp = data.frame(
  var = ordered_component,
  Poisson = varImpStandard(summary(best.poisReg)$coefficients[-1, 3], ordered_component),
  NB = varImpStandard(summary(best.negBino)$coefficients[-1, 3], ordered_component),
  glmnet = varImpStandard(best.glmnet.imp ,ordered_component),
  BRT = varImpStandard(best.gbm.inf, ordered_component),
  stringsAsFactors = FALSE
)

thres = .5  # only include if rel.inf > thres% for any model
varImp_trunc = varImp[apply(varImp[,-1],1,function(x) any(x>thres)),]


library(tidyr)
VI = tidyr::gather(varImp_trunc,model,importance,-var)  # convert to long format


#- order variables by ordered_componenet and replace 0 with NA
VI$var= factor(VI$var,levels=rev(ordered_component))
VI$model = factor(VI$model,levels=plot.info$model)
VI[VI==0] = NA

ggplot(VI, aes(y=var, x=importance)) + 
  geom_point(aes(shape = model, col = model,size=model)) + 
  dot_theme +
  ylab("") + xlab("variable importance") + 
  scale_fill_manual(values=plot.info$col) + 
  scale_shape_manual(values=plot.info$pch) + 
  scale_size_manual(values=3*plot.info$pt.cex) +
  theme(legend.key=element_blank())
  # geom_vline(xintercept = 0)
  #theme(legend.position = c(1.01, -.01), legend.justification=c(1,0), 


ggsave(filename=file.path(plotDir,"Fig3.pdf"),width=7.5,height=7.5)


arrange(varImp_trunc,desc(BRT))

varImp_trunc %>% arrange(desc(BRT)) %>% mutate_each(funs(round),Poisson,NB,glmnet,BRT)
tt=0
varImp_trunc %>% mutate(total=  (NB>tt) + (glmnet>tt) + (BRT>tt)) %>%   # (Poisson>tt)
  arrange(desc(total)) %>% mutate_each(funs(round),Poisson,NB,glmnet,BRT)
##############################################################




#--------------------------------------------
# Figure 4. Interaction effect
#--------------------------------------------
#optimal.gbm = gbm(fmla.offset, data = a, distribution = "poisson", n.trees = 3700, interaction.depth = 4)
#vars = attr(terms(model.frame(fmla.offset, a)),"term.labels")
#p = length(vars)
p = length(ordered_component)

interaction.mag = matrix(0, p, p)

for(i in 1:(p-1)) {
  for(j in (i+1):p) {
    #cat(i, j, "\n")
    var.i = ordered_component[i]
    var.j = ordered_component[j]
    interaction.mag[i, j] = interact.gbm(best.gbm, data = a, i.var = c(var.i, var.j))
  }
}

colnames(interaction.mag) = ordered_component
rownames(interaction.mag) = ordered_component

interaction.mag[is.na(interaction.mag)] = 0


#-- use corrplot()
# order by BRT variable importance
library(dplyr)
imp.vars = arrange(varImp_trunc,desc(BRT)) %>% dplyr::select(var) %>% .[,1] %>% as.character
#imp.vars = as.character(varImp_trunc$var)  # use important variables
i2 = interaction.mag[imp.vars,imp.vars]
i2 = i2 + t(i2)

pdf(file=file.path(plotDir,"Fig4.pdf"),width=5,height=5)
  library(corrplot)
  #par(oma=c(.5,0,4.5,1),xpd=NA)
  #corrplot(i2,is.corr=FALSE,type="upper",tl.col="black",tl.cex=.8,
  #         cl.cex=.7,cl.align.text="l",ylim=c(0,.5),diag=TRUE)
  par(oma=c(.5,0,.5,0),xpd=NA)
  corrplot(i2,is.corr=FALSE,type="full",tl.col="black",tl.cex=.8,
           cl.cex=.7,cl.align.text="l",ylim=c(0,.5),diag=TRUE)
  
  #- Add variable importance score to diagonal
  vi = round(varImp_trunc$BRT[match(imp.vars,varImp_trunc$var)])
  text(1:length(vi),length(vi):1,labels=vi,col="brown",cex=.7)  

dev.off()
  




i2 = interaction.mag[imp.vars,imp.vars]
ij = which(i2 >= 0.10, arr.ind=TRUE)
arrange(data.frame(ij,score=i2[ij]),desc(score)) %>% 
  mutate(i=rownames(i2)[row],j=colnames(i2)[col])

sort(colSums(i2),decreasing=TRUE)
sort(apply(i2,2,function(x) sum(x[x>.10])),decreasing=TRUE)



# BRT.vars = c("Terrain", "SkewAngle", "PedCross", "NextPIDist", "MedWidth", 
#              "Long", "LegWidth", "LegRtType", "Lat", "IntCat") # from the BRT variable importances
# 
# BRT.vars = c("Terrain", "SkewAngle", "PedCross", "NextPIDist", "MedWidth", "LegWidth", "LegRtType", "IntCat")
# interaction.index = which(vars %in% BRT.vars)
# BRT.vars = vars[interaction.index]
# 
# my_palette <- colorRampPalette(c("white", "blue"))(n = 299)
# 
# #par(cex.main=0.5)
# 
# # save as portrait 12 by 6
# heatmap.2(interaction.mag[BRT.vars, BRT.vars], dendrogram = "none", Colv = FALSE,
#           Rowv = FALSE, key = T, keysize = 1, density.info="none", trace="none",
#           key.xlab = "Scale of H-Statistic", key.title = NA,
#           cexRow = 1.4, cexCol = 1.4, margins = c(7, 9), col = my_palette,
#           srtCol=40,revC=TRUE)
#           #main = "Two-way factor interactions")
# 
# 



#--------------------------------------------
# Figure 5. SkewAngle Interactions
#--------------------------------------------

  data = model.frame(fmla.offset,a)  # make correct transformations from fmla
  offset = as.vector(model.offset(data))
  y = model.response(data)    # response vector
  vars = attr(terms(data),"term.labels") # predictor variables
  p = length(vars) 
  
library(mgcv)
  ff = gam(y~te(SkewAngle,LegWidth),family="poisson",offset=offset,data=a)
  plot(ff,scheme=0,se=FALSE,main='')
  plot(ff,scheme=2,se=FALSE,main='')  
  plot(ff,scheme=1,main="crash rate")

  
  




