## preprocessing

cp.poisReg = 1:max.complex
cp.negBino = 1:max.complex
# cp.glmnet = cp.rescale(-(attr(mu.cv.glmnet, "lambda")), 1, max.complex) # note this is reversed lambda
cp.glmnet = attr(mu.cv.glmnet.offset, "lambda")
cp.glmnet = -log(cp.glmnet/max(cp.glmnet))
cp.glmnet = cp.glmnet + (10-max(cp.glmnet))
cp.gbm = tree.seq/1000
cp.tree = cp.seq*200



par(mfrow = c(3, 1), mar=c(6,5,3,4), oma=c(0,2,2,0))

plot(1:max.complex, mae(mu.cv.poisReg.offset, Y), type = 'b', pch = 15, cex = 1.2, 
     ylim=c(1.15, 1.51), col=6, lwd=2, ylab="MAE",  xlab="model complexity",
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.2)
lines(1:max.complex, mae(mu.cv.negBino.offset, Y), type = 'b', col=2, pch = 16, cex = 1.2, lwd=2)
lines(cp.glmnet, mae(mu.cv.glmnet.offset, Y),col=3, lwd=2)
lines(cp.tree, mae(mu.cv.tree.offset, Y), col='orange',lty=4, lwd=2)
lines(cp.gbm,mae(gbm_4.offset,Y),col=4,lty=2, lwd=2)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "Trees", "BRT"), 
       col = c(6, 2, 3, 'orange', 4), pch=c(15,16,NA_integer_,NA_integer_,NA_integer_),
       lty=c(1,1,1,4,2),lwd = 2, text.font = 3, ncol = 5, pt.cex=c(1.2,1.2,1,1,1), cex = 1.2)
title("(a) Prediction error measured by MAE", cex.main = 1.2)


plot(1:max.complex, mse(mu.cv.poisReg.offset, Y), type = 'b', pch = 15, cex = 1.2, 
     col=6,lwd=2, ylab="MSE", ylim=c(5.7, 12), xlab="model complexity",
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.2)
lines(1:max.complex, mse(mu.cv.negBino.offset, Y), type = 'b', col=2, pch = 16, cex = 1.2,lwd=2)
lines(cp.glmnet, mse(mu.cv.glmnet.offset, Y),col=3,lwd=2)
lines(cp.tree, mse(mu.cv.tree.offset, Y), col='orange',lty=4,lwd=2)
lines(cp.gbm,mse(gbm_4.offset,Y),col=4,lty=2,lwd=2)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "Trees", "BRT"), 
       col = c(6, 2, 3, 'orange', 4), pch=c(15,16,NA_integer_,NA_integer_,NA_integer_),
       lty=c(1,1,1,4,2),lwd = 2, text.font = 3, ncol = 5, pt.cex=c(1.2,1.2,1,1,1), cex = 1.2)
title("(b) Prediction error measured by MSE", cex.main = 1.2)



plot(1:max.complex, mlogL(mu.cv.poisReg.offset, Y), type = 'b', pch = 15, cex = 1.2, 
     col=6,lwd=2, ylab="-mlogL", ylim = c(1.57, 2.2),  xlab="model complexity",
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.2)
lines(1:max.complex, mlogL(mu.cv.negBino.offset, Y), type = 'b', pch = 16, cex = 1.2, col=2,lwd=2)
lines(cp.glmnet, mlogL(mu.cv.glmnet.offset, Y),col=3,lwd=2)
lines(cp.tree, mlogL(mu.cv.tree.offset, Y), col='orange',lty=4,lwd=2)
lines(cp.gbm,mlogL(gbm_4.offset,Y),col=4,lty=2,lwd=2)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "Trees", "BRT"), 
       col = c(6, 2, 3, 'orange', 4), pch=c(15,16,NA_integer_,NA_integer_,NA_integer_),
       lty=c(1,1,1,4,2),lwd = 2, text.font = 3, ncol = 5, pt.cex=c(1.2,1.2,1,1,1), cex = 1.2)
#title("(c) Prediction error measured by -mlogL", cex.main = 1.2)
title("Prediction error measured by -mlogL", cex.main = 1.2)

## output both transformed and original optimal complexities
# original scale
(1:max.complex)[which.min( mlogL(mu.cv.poisReg.offset , Y))]
(1:max.complex)[which.min( mlogL(mu.cv.negBino.offset , Y))]
(attr(mu.cv.glmnet.offset, "lambda"))[which.min(mlogL(mu.cv.glmnet.offset , Y))]
tree.seq[which.min(mlogL(gbm_4.offset , Y))]
cp.seq[which.min(mlogL(mu.cv.tree.offset , Y))]

# original scale
(1:max.complex)[which.min( mlogL(mu.cv.poisReg.offset , Y))]
(1:max.complex)[which.min( mlogL(mu.cv.negBino.offset , Y))]
cp.glmnet[which.min(mlogL(mu.cv.glmnet.offset , Y))]
cp.gbm[which.min(mlogL(gbm_4.offset , Y))]
cp.tree[which.min(mlogL(mu.cv.tree.offset , Y))]



