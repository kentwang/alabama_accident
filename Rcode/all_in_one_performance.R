## preprocessing

cp.poisReg = 1:max.complex
cp.negBino = 1:max.complex
cp.glmnet = cp.rescale(-(attr(mu.cv.glmnet, "lambda")), 1, max.complex) # note this is reversed lambda
cp.gbm = cp.rescale(tree.seq, 1, max.complex)
cp.tree = cp.rescale(cp.seq, 1, max.complex)



par(mfrow = c(3, 1))

plot(1:max.complex, mae(mu.cv.poisReg.offset, Y), type = 'l', ylim=c(1.15, 1.6), col=6, ylab="MAE",  xlab="cp*", main="Comparison of families on Mean Absolute Error")
lines(1:max.complex, mae(mu.cv.negBino.offset, Y), col=2)
lines(cp.glmnet, mae(mu.cv.glmnet.offset, Y),col=3)
lines(cp.gbm,mae(gbm_4.offset,Y),col=4)
lines(cp.tree, mae(mu.cv.tree.offset, Y), col=5)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "BRT", "Trees"), 
       col = c(6, 2, 3, 4, 5), lwd = 2, text.font = 3, ncol = 2)


plot(1:max.complex, mse(mu.cv.poisReg.offset, Y), type = 'l', col=6, ylab="MSE", ylim=c(5.7, 9), xlab="cp*", main="Comparison of families on Mean Square Error")
lines(1:max.complex, mse(mu.cv.negBino.offset, Y), col=2)
lines(cp.glmnet, mse(mu.cv.glmnet.offset, Y),col=3)
lines(cp.gbm,mse(gbm_4.offset,Y),col=4)
lines(cp.tree, mse(mu.cv.tree.offset, Y), col=5)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "BRT", "Trees"), 
       col = c(6, 2, 3, 4, 5), lwd = 2, text.font = 3, ncol = 2)

plot(1:max.complex, mlogL(mu.cv.poisReg.offset, Y), type = 'l', col=6, ylab="-mlogL", ylim = c(1.57, 2.0),  xlab="cp*", main="Comparison of families on Negative Maximumzed Log-likelihood")
lines(1:max.complex, mlogL(mu.cv.negBino.offset, Y), col=2)
lines(cp.glmnet, mlogL(mu.cv.glmnet.offset, Y),col=3)
lines(cp.gbm,mlogL(gbm_4.offset,Y),col=4)
lines(cp.tree, mlogL(mu.cv.tree.offset, Y), col=5)
legend("topleft" ,legend = c("Poisson", "NB", "glmnet", "BRT", "Trees"), 
       col = c(6, 2, 3, 4, 5), lwd = 2, text.font = 3, ncol = 2)