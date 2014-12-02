## preprocessing

cp.poisReg = 1:max.complex
cp.negBino = 1:max.complex
# cp.glmnet = cp.rescale(-(attr(mu.cv.glmnet, "lambda")), 1, max.complex) # note this is reversed lambda
cp.glmnet = attr(mu.cv.glmnet, "lambda")
cp.glmnet = -log(cp.glmnet/max(cp.glmnet))
cp.glmnet = cp.glmnet + (10-max(cp.glmnet))
cp.gbm = tree.seq/1000
cp.tree = cp.seq*200



par(mfrow = c(3, 1), mar=c(5,4,1,2), oma=c(0,1,2,0))

plot(1:max.complex, mae(mu.cv.poisReg.offset, Y), type = 'b', pch = 15, cex = 1.2, ylim=c(1.15, 1.51), col=6, lwd=2, ylab="MAE",  xlab="mcp*")
lines(1:max.complex, mae(mu.cv.negBino.offset, Y), type = 'b', col=2, pch = 16, cex = 1.2, lwd=2)
lines(cp.glmnet, mae(mu.cv.glmnet.offset, Y),col=3, lwd=2)
lines(cp.tree, mae(mu.cv.tree.offset, Y), col='orange',lty=4, lwd=2)
lines(cp.gbm,mae(gbm_4.offset,Y),col=4,lty=2, lwd=2)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "Trees", "BRT"), 
       col = c(6, 2, 3, 'orange', 4), pch=c(15,16,NA_integer_,NA_integer_,NA_integer_),
       lty=c(1,1,1,4,2),lwd = 2, text.font = 3, ncol = 5, pt.cex=c(1.2,1.2,1,1,1))


plot(1:max.complex, mse(mu.cv.poisReg.offset, Y), type = 'b', pch = 15, cex = 1.2, col=6,lwd=2, ylab="MSE", ylim=c(5.7, 12), xlab="mcp*")
lines(1:max.complex, mse(mu.cv.negBino.offset, Y), type = 'b', col=2, pch = 16, cex = 1.2,lwd=2)
lines(cp.glmnet, mse(mu.cv.glmnet.offset, Y),col=3,lwd=2)
lines(cp.tree, mse(mu.cv.tree.offset, Y), col='orange',lty=4,lwd=2)
lines(cp.gbm,mse(gbm_4.offset,Y),col=4,lty=2,lwd=2)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "Trees", "BRT"), 
       col = c(6, 2, 3, 'orange', 4), pch=c(15,16,NA_integer_,NA_integer_,NA_integer_),
       lty=c(1,1,1,4,2),lwd = 2, text.font = 3, ncol = 5, pt.cex=c(1.2,1.2,1,1,1))


plot(1:max.complex, mlogL(mu.cv.poisReg.offset, Y), type = 'b', pch = 15, cex = 1.2, col=6,lwd=2, ylab="-mlogL", ylim = c(1.57, 2.2),  xlab="mcp*")
lines(1:max.complex, mlogL(mu.cv.negBino.offset, Y), type = 'b', pch = 16, cex = 1.2, col=2,lwd=2)
lines(cp.glmnet, mlogL(mu.cv.glmnet.offset, Y),col=3,lwd=2)
lines(cp.tree, mlogL(mu.cv.tree.offset, Y), col='orange',lty=4,lwd=2)
lines(cp.gbm,mlogL(gbm_4.offset,Y),col=4,lty=2,lwd=2)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "Trees", "BRT"), 
       col = c(6, 2, 3, 'orange', 4), pch=c(15,16,NA_integer_,NA_integer_,NA_integer_),
       lty=c(1,1,1,4,2),lwd = 2, text.font = 3, ncol = 5, pt.cex=c(1.2,1.2,1,1,1))


title("Comparison of families on MAE, MSE, and -mlogL", outer=TRUE)