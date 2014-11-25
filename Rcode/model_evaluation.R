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

mu.cv.tree = cv.treePois.old(fmla,data=a,fold=fold)
mu.cv.tree.offset = cv.treePois.old(fmla.offset,data=a,fold=fold)

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
#-- Use cross-validation (the cross validation framework has been corrected)
#   Compare poisReg, negBino, gbmPois, glmnet, tree
################################################################################
fold = cvfolds(nrow(a),k=20,seed=9122014)  # get cv partition info fold 20 and 5

Y = a$X5YrCrashCount  # response


## modified cross validation for Poisson Regression
# Todo: NA should be removed
max.complex = 10 # try 15?
mu.cv.poisReg <- cv.poisReg(fmla, data = a, fold = fold, max.complex = max.complex)
mu.cv.poisReg.offset <- cv.poisReg(fmla.offset, data = a, fold= fold, max.complex = max.complex)

par(mfrow = c(3, 1))
plot(1:max.complex, mae(mu.cv.poisReg , Y), typ='l', col=3, ylab="MAE", xlab="# of predictors")
lines(1:max.complex, mae(mu.cv.poisReg.offset, Y), lty=2, col=3)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(3, 3), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Poisson regression using MAE")

# regression tree performance using MSE
plot(1:max.complex, mse(mu.cv.poisReg , Y), typ='l', col=4, ylab="MSE", xlab="# of predictors")
lines(1:max.complex, mse(mu.cv.poisReg.offset, Y), lty=2, col=4)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(4, 4), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Poisson regression using MAE")

# regression tree performance using mlogL
plot(1:max.complex, mlogL(mu.cv.poisReg, Y), typ='l', col=6, ylab="mlogL",  xlab="# of predictors")
lines(1:max.complex, mlogL(mu.cv.poisReg.offset, Y), lty=2, col=6)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(6, 6), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Poisson regression using MAE")


mu.cv.poisReg.old <- suppressWarnings(cv.poisReg.old(fmla, data = a, fold = fold))
mu.cv.poisReg.offset.old <- suppressWarnings(cv.poisReg.old(fmla.offset, data = a, fold = fold))



## modified cross validation for Negative Binomial regression
max.complex = 10
mu.cv.negBino <- cv.negBino(fmla, data = a, fold = fold, max.complex = max.complex)
mu.cv.negBino.offset <- cv.negBino(fmla.offset, data = a, fold= fold, max.complex = max.complex)

par(mfrow = c(3, 1))
# regression tree performance using MAE
plot(1:max.complex, mae(mu.cv.negBino , Y), typ='l', col=3, ylab="MAE", xlab="# of predictors")
lines(1:max.complex, mae(mu.cv.negBino.offset, Y), lty=2, col=3)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(3, 3), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Negative Binomial regression using MAE")

# regression tree performance using MSE
plot(1:max.complex, mse(mu.cv.negBino , Y), typ='l', col=4, ylab="MSE", xlab="# of predictors")
lines(1:max.complex, mse(mu.cv.negBino.offset, Y), lty=2, col=4)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(4, 4), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Negative Binomial regression using MSE")

# regression tree performance using mlogL
plot(1:max.complex, mlogL(mu.cv.negBino, Y), typ='l', col=6, ylab="mlogL",  xlab="# of predictors")
lines(1:max.complex, mlogL(mu.cv.negBino.offset, Y), lty=2, col=6)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(6, 6), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Negative Binomial regression using mlogL")


mu.cv.negBino <- suppressWarnings(cv.negBino.old(fmla, data = a, fold = fold))
mu.cv.negBino.offset <- suppressWarnings(cv.negBi
                                         no.old(fmla.offset, data = a, fold = fold))


## modified cross validation for glmnet
mu.cv.glmnet = cv.glmnetPois(fmla,data=a,fold=fold)
mu.cv.glmnet.offset = cv.glmnetPois(fmla.offset,data=a,fold=fold)

par(mfrow = c(3, 1))
# regression tree performance using MAE
plot(attr(mu.cv.glmnet, "lambda"), mae(mu.cv.glmnet , Y), typ='l', xlim=c(0, 1.0), ylim=c(1.29, 1.5), col=3, ylab="MAE", xlab=expression(~lambda))
lines(attr(mu.cv.glmnet.offset, "lambda"), mae(mu.cv.glmnet.offset, Y), lty=2, col=3)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(3, 3), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Glmnet Poisson using MAE")

# regression tree performance using MSE
plot(attr(mu.cv.glmnet, "lambda"), mse(mu.cv.glmnet , Y), typ='l', xlim=c(0, 1.0), ylim =c(6.95, 10), col=4, ylab="MSE", xlab=expression(~lambda))
lines(attr(mu.cv.glmnet.offset, "lambda"), mse(mu.cv.glmnet.offset, Y), lty=2, col=4)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(4, 4), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Glmnet Poisson using MSE")

# regression tree performance using mlogL
plot(attr(mu.cv.glmnet, "lambda"), mlogL(mu.cv.glmnet , Y), typ='l', xlim=c(0, 1.0), ylim =c(1.73, 2.1), col=6, ylab="mlogL",  xlab=expression(~lambda))
lines(attr(mu.cv.glmnet.offset, "lambda"), mlogL(mu.cv.glmnet.offset, Y), lty=2, col=6)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(6, 6), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of Glmnet Poisson using mlogL")

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


par(mfrow = c(3, 1))
# Plot gbm tree performance in MAE
plot(tree.seq,mae(gbm_3,Y),typ='o',ylim=c(1.18,1.5),col=3,ylab="MAE")
lines(tree.seq,mae(gbm_3.offset,Y),lty=3,col=3)
lines(tree.seq,mae(gbm_1,Y),col=1)
lines(tree.seq,mae(gbm_1.offset,Y),lty=3,col=1)
lines(tree.seq,mae(gbm_2,Y),col=2)
lines(tree.seq,mae(gbm_2.offset,Y),lty=3,col=2)
lines(tree.seq,mae(gbm_4,Y),col=4)
lines(tree.seq,mae(gbm_4.offset,Y),lty=3,col=4)

abline(h=mae(gbm.old_3,Y),col=5)
abline(h=mae(gbm.old_3.offset,Y),col=5,lty=3)
title("GBM Comparison Using Mean Absolute Error")
legend("topleft" ,legend = c("Dept 3", "Dept 1", "Dept 2", "Dept 4",
                             "Dept 3 offset", "Dept 1 offset", "Dept 2 offset", "Dept 4 offset"),
       lty = c(1, 1, 1, 1, 3, 3, 3, 3), col = c(3, 1, 2, 4, 3, 1, 2, 4), lwd = 2,
       ncol = 2, text.font = 3, cex = 0.8)



# Note: Repeat using mean squared error, likelihood, etc. 
# plot gbm tree performance in MSE
plot(tree.seq,mse(gbm_3,Y),type="l",ylim=c(5.8, 8.5), col=3,ylab="MSE")
lines(tree.seq,mse(gbm_3.offset,Y),lty=3,col=3)
lines(tree.seq,mse(gbm_1,Y),col=1)
lines(tree.seq,mse(gbm_1.offset,Y),lty=3,col=1)
lines(tree.seq,mse(gbm_2,Y),col=2)
lines(tree.seq,mse(gbm_2.offset,Y),lty=3,col=2)
lines(tree.seq,mse(gbm_4,Y),col=4)
lines(tree.seq,mse(gbm_4.offset,Y),lty=3,col=4)

abline(h=mse(gbm.old_3,Y),col=5)
abline(h=mse(gbm.old_3.offset,Y),col=5,lty=3)
title("GBM Comparison Using Mean Square Error")
legend("topleft" ,legend = c("Dept 3", "Dept 1", "Dept 2", "Dept 4",
                             "Dept 3 offset", "Dept 1 offset", "Dept 2 offset", "Dept 4 offset"),
       lty = c(1, 1, 1, 1, 3, 3, 3, 3), col = c(3, 1, 2, 4, 3, 1, 2, 4), lwd = 2,
       ncol = 2, text.font = 3, cex = 0.8)

# plot gbm tree performance in likelihood
plot(tree.seq,mlogL(gbm_3,Y),type="l", ylim = c(1.3, 1.82), col=3,ylab="mlogL")
lines(tree.seq,mlogL(gbm_3.offset,Y),lty=3,col=3)
lines(tree.seq,mlogL(gbm_1,Y),col=1)
lines(tree.seq,mlogL(gbm_1.offset,Y),lty=3,col=1)
lines(tree.seq,mlogL(gbm_2,Y),col=2)
lines(tree.seq,mlogL(gbm_2.offset,Y),lty=3,col=2)
lines(tree.seq,mlogL(gbm_4,Y),col=4)
lines(tree.seq,mlogL(gbm_4.offset,Y),lty=3,col=4)

abline(h=mlogL(gbm.old_3,Y),col=5)
abline(h=mlogL(gbm.old_3.offset,Y),col=5,lty=3)
title("GBM Comparison Using Mean Loglikelihood")
legend("topleft" ,legend = c("Dept 3", "Dept 1", "Dept 2", "Dept 4",
                             "Dept 3 offset", "Dept 1 offset", "Dept 2 offset", "Dept 4 offset"),
       lty = c(1, 1, 1, 1, 3, 3, 3, 3), col = c(3, 1, 2, 4, 3, 1, 2, 4), lwd = 2,
       ncol = 2, text.font = 3, cex = 0.8)


#### modified cross validation of trees
# set.seed(11102014)
# cp.seq = sort(rgamma(100, 1, 20), decreasing = T) # gamma distribution seems to be what I want for the seach of cp (right skewed)
cp.seq = seq(0, 0.05, by = 0.0001) # need tuning the parameter
# hist(cp.seq, breaks = 20)
# plot(cp.seq)
length(cp.seq)

mu.cv.tree = cv.treePois(fmla, data = a, fold = fold, cp.seq = cp.seq, show.pb=TRUE)
mu.cv.tree.offset = cv.treePois(fmla.offset, data = a, fold = fold, cp.seq = cp.seq, show.pb=TRUE)

par(mfrow = c(3, 1))

# regression tree performance using MAE
plot(cp.seq, mae(mu.cv.tree, Y), typ='l', col=3, ylim = c(1.28, 1.5), ylab="MAE")
lines(cp.seq, mae(mu.cv.tree.offset, Y), lty=2, col=3)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(3, 3), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of regression tree using MAE")

# regression tree performance using MSE
plot(cp.seq, mse(mu.cv.tree, Y), typ='l', col=4, ylim = c(6.5, 8.5), ylab="MSE")
lines(cp.seq, mse(mu.cv.tree.offset, Y), lty=2, col=4)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(4, 4), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of regression tree using MSE")

# regression tree performance using mlogL
plot(cp.seq, mlogL(mu.cv.tree, Y), typ='l', col=6, ylim = c(1.7, 2.0), ylab="mlogL")
lines(cp.seq, mlogL(mu.cv.tree.offset, Y), lty=2, col=6)
legend("topleft" ,legend = c("No offset", "Offset"), lty = c(1, 3), col = c(6, 6), lwd = 2, text.font = 3, cex = 0.8)
title("Performace of regression tree using mlogL")


mu.cv.tree.old = cv.treePois.old(fmla,data=a,fold=fold)
mu.cv.tree.offset.old = cv.treePois.old(fmla.offset,data=a,fold=fold)


##### following no need to run
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
