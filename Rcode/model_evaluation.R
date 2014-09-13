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


#-- Evaluation (mean absolut error)
Y = a$X5YrCrashCount  # response

score = data.frame(tree=mae(m.tree,Y),glm=mae(m.glm,Y),step=mae(m.step,Y),
                    glmnet=min(mae(m.glmnet,Y)))

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

