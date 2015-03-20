# This program fits the best model for each family using full data
## poisson
max.complex = 10
mu.cv.poisReg.offset <- cv.poisReg(fmla.offset, data = a, fold= fold, max.complex = max.complex)
plot(1:max.complex, mlogL(mu.cv.poisReg.offset , Y), typ='l', col=3, ylab="mlogL", xlab="# of predictors")


## nagative binomial
mu.cv.negBino.offset <- cv.negBino(fmla.offset, data = a, fold= fold, max.complex = max.complex)
plot(1:max.complex, mlogL(mu.cv.negBino.offset , Y), typ='l', col=3, ylab="mlogL", xlab="# of predictors")


## glmnet poisson
mu.cv.glmnet.offset = cv.glmnetPois(fmla.offset,data=a,fold=fold)
plot(attr(mu.cv.glmnet.offset, "lambda"), mlogL(mu.cv.glmnet.offset , Y), typ='l', xlim=c(0, 1.0), col=3, ylab="mlogL", xlab=expression(~lambda))

## gbm trees
tree.seq = seq(500,10000,by=100)
gbm_4.offset <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                           tree.seq = tree.seq, interaction.depth=4,
                           shrinkage = .005)
plot(tree.seq,mlogL(gbm_4.offset,Y),typ='l',col=3,ylab="mlogL")


## regression trees
cp.seq = seq(0, 0.05, by = 0.0001) # need tuning the parameter
mu.cv.tree.offset = cv.treePois(fmla.offset, data = a, fold = fold, cp.seq = cp.seq, show.pb=TRUE)
plot(cp.seq, mlogL(mu.cv.tree.offset, Y), typ='l', col=3, ylab="mlogL")


## optimal complexity parameters
cp.poisReg = (1:max.complex)[which.min( mlogL(mu.cv.poisReg.offset , Y))]
cp.negBino = (1:max.complex)[which.min( mlogL(mu.cv.negBino.offset , Y))]
cp.glmnet = (attr(mu.cv.glmnet.offset, "lambda"))[which.min(mlogL(mu.cv.glmnet.offset , Y))]
cp.gbm = tree.seq[which.min(mlogL(gbm_4.offset , Y))]
cp.tree = cp.seq[which.min(mlogL(mu.cv.tree.offset , Y))]

## fit the best model out of each family
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

# This is a showcase of regress trees
tree.test = rpart(fmla.offset, data = a, method = "poisson", cp = 0.02,xval=0, minbucket=3)
plot(tree.test, uniform=TRUE, main="Regression Tree for Accident Count")
text(tree.test, cex = 1.5)

## importance comparison
fmla.string = strsplit(as.character(fmla.offset)[3], " \\+ ")[[1]]
fmla.string = fmla.string[fmla.string != "offset(log(Traffic))"]
allVariables = fmla.string[order(fmla.string)]
allVariables = allVariables[-c(11, 12)] # remove log(TotalAADT) and log(TotalT5YearInM

# processing of glmnet
# Todo standardization of coefficient by sd(X)
best.glmnet.imp = as.vector(best.glmnet$beta)/apply(X, 2, sd) # standardized by the coded model matrix
names(best.glmnet.imp) = rownames(best.glmnet$beta)

# processing of gbm
best.gbm.inf = best.gbm$rel.inf
names(best.gbm.inf) = best.gbm$var

# reorder the importance using the order of the component graph
ordered_component = c("IntCat", "LegRtType", "AreaType", "LegWidth", "OneWay", "Lat",
                      "Terrain", "Long", "NextPIDist", "Lighting", "LegTCType", "LegSpeed",
                      "MedWidth", "LegType", "OffsetDist", "RTMoveCtrl", "Offset",
                      "LTWidth", "IntTCType", "LTOffset", "NumLanes", "SightRt",
                      "LTLanes", "LTLnLength", "PedCross", "RTWidth", "SightLt",
                      "PaveType", "TurnProhib", "MedType", "RTLanes", "RTChannel",
                      "MergeLanes", "NumberLegs", "Rumble", "SkewAngle", "RTLnLength")

comp_order = match(ordered_component, allVariables) # indices from allvariabels to ordered_component


## combind importances
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



ggplot(importances, aes(factor(var), importance, fill = model)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.4), width=.8) +
  theme(text = element_text(size = 16)) + 
  coord_flip() + 
  scale_x_discrete(name="Intersection Factors") +
  scale_fill_discrete(breaks=model.names)

# use cleveland dot plot
# dotchart(VADeaths, main = "Death Rates in Virginia - 1940")

# cleveland dot chart group by model
importances_cdp = t(matrix(varImp, ncol = 37, byrow = TRUE))
colnames(importances_cdp) = model.names
rownames(importances_cdp) = ordered_component

dotchart(importances_cdp, cex = 0.5)

# cleveland dot chart group by factors
importances_cdp = matrix(varImp, ncol = 37, byrow = TRUE)
rownames(importances_cdp) = model.names
colnames(importances_cdp) = ordered_component

dotchart(importances_cdp, cex = 0.5)

# <<<<<<< HEAD
# =======
### Checking interactions using interaction.gbm
# todo: all two way interactions top 10
optimal.gbm = gbm(fmla.offset, data = a, distribution = "poisson", n.trees = 3700, interaction.depth = 4)
vars = attr(terms(model.frame(fmla.offset, a)),"term.labels")
p = length(vars)

interaction.mag = matrix(0, p, p)

for(i in 1:(p-1)) {
  for(j in (i+1):p) {
    cat(i, j, "\n")
    interaction.mag[i, j] = interact.gbm(optimal.gbm, data = a, i.var = c(i, j), n.trees = 3700)
    # cat(vars[i], vars[j], interact.gbm(optimal.gbm, data = a, i.var = c(i, j), n.trees = 3700), "\n")
  }
}

colnames(interaction.mag) = vars
rownames(interaction.mag) = vars

interaction.mag[is.na(interaction.mag)] = 0
#interaction.mag[interaction.mag == 0] = 10^-6

#graph too big. pick some representative variables
#heatmap(interaction.mag, Colv = NA, Rowv= NA)
#interaction.index = sort(c(37, 36, 33, 29, 24, 15, 8, 4))
#heatmap(interaction.mag[interaction.index, interaction.index])

BRT.vars = c("Terrain", "SkewAngle", "PedCross", "NextPIDist", "MedWidth", 
             "Long", "LegWidth", "LegRtType", "Lat", "IntCat") # from the BRT variable importances

BRT.vars = c("Terrain", "SkewAngle", "PedCross", "NextPIDist", "MedWidth", "LegWidth", "LegRtType", "IntCat")
interaction.index = which(vars %in% BRT.vars)
BRT.vars = vars[interaction.index]
# heatmap(interaction.mag[BRT.vars, BRT.vars], Colv = NA, Rowv = NA, cexRow = 1, cexCol = 1)

my_palette <- colorRampPalette(c("yellow", "red"))(n = 299)

par(cex.main=0.5)

heatmap.2(interaction.mag[BRT.vars, BRT.vars], dendrogram = "none",
          key = T, keysize = 1.2, density.info="none", trace="none",
          key.xlab = "Scale of H-Statistic", key.title = NA,
          cexRow = 1.2, cexCol = 1.2, margins = c(7, 9), col = my_palette,
          srtCol=45, main = "Two-way factor interactions")
# >>>>>>> 981d79ab312916467e596739eb02f819efef94c8


