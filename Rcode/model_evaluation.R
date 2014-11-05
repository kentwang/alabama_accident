#-------------------------------------------------------------------------------
#- Compare models yielded from model_selection.R
#
#- Todo: -regular cross-validation 1) abs erro 2) likelihood? trees
#-------------------------------------------------------------------------------

#- Some functions
#- Todo: -validating functions

cvDataSplit <- function(data, k) {
  dataList <- list()
  n <- dim(data)[1]
  subn <- ceiling(n / k)
  rid <- sample(1:n)
  for(i in 1:(k - 1)) { # not including the last set
    dataList <- c(dataList, list(data[1:subn, ]))
    rid <- rid[-(1:subn)]
  }
  dataList <- c(dataList, list(data[rid, ]))
  return(dataList)
}







source("RCode/functions.R")  # load required functions
source('Rcode/glm_my.R')

################################################################################
#-- Fit models to full data. Didn't put in offset models yet. Basically, this
#   is a duplicate of what's in model_selection.R
#
#   Note: will be prone to overfitting (especially glm)
################################################################################
m.tree = predict(rpart(fmla,data=a,method="poisson"))

X = model.matrix(fmla,data=a)[,-1]
Y = as.matrix(a[,as.character(fmla[[2]])])
offset = model.offset(model.frame(fmla,data=a))
fit = glmnet(X,Y,offset=offset,family="poisson",alpha=0.8)
m.glmnet = predict(fit,type="response",newx=X,offset=offset,s=fit$lambda)


p1 = glm(fmla, data = a, family = "poisson")
m.glm = predict(p1,type="response")

poisReg <- step(glm(fmla, data = a, family = "poisson"))
m.step = predict(poisReg,type="response")

gbmPois <- gbm(fmla, data = a, distribution = "poisson", shrinkage = 1.001, n.trees = 1000) # shrinkage needs to be tuned
best.iter <- gbm.perf(gbmPois, plot.it = FALSE)
m.gbm <- predict(gbmPois, a, n.trees = best.iter, type = "response")


#-- Evaluation (mean absolut error)
Y = a$X5YrCrashCount  # response

score = data.frame(tree=mae(m.tree,Y),glm=mae(m.glm,Y),step=mae(m.step,Y),
                    glmnet=min(mae(m.glmnet,Y)), gbm = mae(m.gbm, Y))

print(round(score,4))


################################################################################
#-- Use cross-validation
#   Note: rpart uses cv for pruning (or stopping). It may be more reasonable to
#    let tree grow to n number of leaves, then determing this value from our cv
################################################################################
fold = cvfolds(nrow(a),k=20,seed=9122014)  # get cv partition info

#mu.cv.step = cv.poisReg(fmla,data=a,fold=fold)
#mu.cv.step.offset = cv.poisReg(fmla.offset,data=a,fold=fold)

mu.cv.glmnet = cv.glmnetPois(fmla,data=a,fold=fold)
mu.cv.glmnet.offset = cv.glmnetPois(fmla.offset,data=a,fold=fold)

mu.cv.tree = cv.treePois(fmla,data=a,fold=fold)
mu.cv.tree.offset = cv.treePois(fmla.offset,data=a,fold=fold)

#-- Evaluation (mean absolut error)
Y = a$X5YrCrashCount  # response

score = data.frame(tree=mae(mu.cv.tree,Y),tree.offset=mae(mu.cv.tree.offset,Y),
                   glmnet=min(mae(mu.cv.glmnet,Y)),glmnet.offset=min(mae(mu.cv.glmnet.offset,Y)))

print(round(score,4))

plot(mae(mu.cv.glmnet,Y),ylim=c(min(score),1.6),typ='l')
lines(mae(mu.cv.glmnet.offset,Y),col=2)
abline(h=mae(mu.cv.tree,Y),col=3)
abline(h=mae(mu.cv.tree.offset,Y),col=4)
legend("topright",c('glmnet','glmnet-offset','tree','tree-offset'),col=1:4,lwd=1)


################################################################################
#-- Use cross-validation
#   Compare poisReg, negBino, gbmPois, glmnet, tree
################################################################################
ptm <- proc.time() # time the comparison procedure

fold = cvfolds(nrow(a),k=20,seed=9122014)  # get cv partition info fold 20 and 5

Y = a$X5YrCrashCount  # response

mu.cv.poisReg <- suppressWarnings(cv.poisReg(fmla, data = a, fold = fold))
mu.cv.poisReg.offset <- suppressWarnings(cv.poisReg(fmla.offset, data = a, fold = fold))

mu.cv.negBino <- suppressWarnings(cv.negBino(fmla, data = a, fold = fold))
mu.cv.negBino.offset <- suppressWarnings(cv.negBi
                                         no(fmla.offset, data = a, fold = fold))



## Boosted Trees
#  - number according to interaction depth. 
#  - Set shrinkage at .005 so converge faster than default of .001



#mu.cv.gbmPois <- suppressWarnings(cv.gbmPois(fmla, data = a, fold = fold))
#mu.cv.gbmPois.offset <- suppressWarnings(cv.gbmPois(fmla.offset, data = a, fold = fold))

tree.seq = seq(500,10000,by=100)
gbm_3 <- cv.gbmPois(fmla, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=3,
                            shrinkage = .005)

gbm_3.offset <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=3,
                            shrinkage = .005)

gbm.old_3 <- cv.gbmPois.old(fmla, data = a, fold = fold,
                                    max.trees=10000, interaction.depth=3,
                                    shrinkage = .005)

gbm.old_3.offset <- cv.gbmPois.old(fmla.offset, data = a, fold = fold,
                                       max.trees=10000,interaction.depth=3,
                                       shrinkage = .005)

gbm_2 <- cv.gbmPois(fmla, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=2,
                            shrinkage = .005)

gbm_2.offset <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=2,
                            shrinkage = .005)

gbm_1 <- cv.gbmPois(fmla, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=1,
                            shrinkage = .005)

gbm_1.offset <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=1,
                            shrinkage = .005)

gbm_4 <- cv.gbmPois(fmla, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=4,
                            shrinkage = .005)

gbm_4.offset <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=4,
                            shrinkage = .005)


# Plot tree performance
plot(tree.seq,mae(gbm_3,Y),typ='o',ylim=c(1.2,1.4),col=3,ylab="MAE")
lines(tree.seq,mae(gbm_3.offset,Y),lty=3,col=3)
lines(tree.seq,mae(gbm_1,Y),col=1)
lines(tree.seq,mae(gbm_1.offset,Y),lty=3,col=1)
lines(tree.seq,mae(gbm_2,Y),col=2)
lines(tree.seq,mae(gbm_2.offset,Y),lty=3,col=2)
lines(tree.seq,mae(gbm_4,Y),col=4)
lines(tree.seq,mae(gbm_4.offset,Y),lty=3,col=4)

abline(h=mae(gbm.old_3,Y),col=5)
abline(h=mae(gbm.old_3.offset,Y),col=5,lty=3)



# Note: Repeat using mean squared error, likelihood, etc. 




# > proc.time() - ptm
# user  system elapsed 
# 1300.86    0.00 1302.19 

mu.cv.glmnet = cv.glmnetPois(fmla,data=a,fold=fold)
mu.cv.glmnet.offset = cv.glmnetPois(fmla.offset,data=a,fold=fold)

mu.cv.tree = cv.treePois(fmla,data=a,fold=fold)
mu.cv.tree.offset = cv.treePois(fmla.offset,data=a,fold=fold)

proc.time() - ptm # calculate the comparison

score = data.frame(poisReg=mae(mu.cv.poisReg,Y),poisReg.offset=mae(mu.cv.poisReg.offset,Y),
                   negBino=mae(mu.cv.negBino,Y),negBino.offset=mae(mu.cv.negBino.offset,Y),
                   gbmPois=mae(mu.cv.gbmPois,Y),gbmPois.offset=mae(mu.cv.gbmPois.offset,Y),
                   tree=mae(mu.cv.tree,Y),tree.offset=mae(mu.cv.tree.offset,Y),
                   glmnet=min(mae(mu.cv.glmnet,Y)),glmnet.offset=min(mae(mu.cv.glmnet.offset,Y)))
print(round(score,4))

plot(1:length(score), score, ylim = c(1.25, 1.55), "h", ylab = "MAE", xlab = "", 
     xaxt = "n", lwd = 2, main = "Comparison of model families on MAE")
points(1:length(score), score, pch = 19)
axis(1, at = 1:10, labels=FALSE)
text(x = 1:10, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
     labels = names(score), srt = 30, adj = 1, xpd =TRUE)
