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
plot(attr(mu.cv.glmnet, "lambda"), mlogL(mu.cv.glmnet.offset , Y), typ='l', xlim=c(0, 1.0), col=3, ylab="mlogL", xlab=expression(~lambda))

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
cp.glmnet = (attr(mu.cv.glmnet, "lambda"))[which.min(mlogL(mu.cv.glmnet.offset , Y))]
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

## importance comparison
fmla.string = strsplit(as.character(fmla.offset)[3], " \\+ ")[[1]]
fmla.string = fmla.string[fmla.string != "offset(log(Traffic))"]
allVariables = fmla.string[order(fmla.string)]

# processing of glmnet
# Todo standardization of coefficient by sd(X)
best.glmnet.imp = as.vector(best.glmnet$beta)/apply(X, 2, sd) # standardized by the coded model matrix
names(best.glmnet.imp) = rownames(best.glmnet$beta)

# processing of gbm
best.gbm.inf = best.gbm$rel.inf
names(best.gbm.inf) = best.gbm$var

## combind importances
varImp <- as.vector(c(varImpStandard(summary(best.poisReg)$coefficients[-1, 3], allVariables),
                      varImpStandard(summary(best.negBino)$coefficients[-1, 3], allVariables),
                    varImpStandard(best.glmnet.imp ,allVariables),
                    varImpStandard(best.tree$variable.importance, allVariables),
                    varImpStandard(best.gbm.inf, allVariables)))

model.names <- c("Poisson", "NB",
                 "glmnet", "Trees", "BRT")

importances <- as.data.frame(cbind(rep(allVariables, length(model.names)),
                                   rep(model.names, each = length(allVariables))))

importances$importance <- varImp
names(importances) <- c("var", "model", "importance")

ggplot(importances, aes(factor(var), importance, fill = model)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(text = element_text(size = 16)) + 
  coord_flip() + 
  scale_x_discrete(name="Intersection Factors") +
  scale_fill_discrete(breaks=model.names)



