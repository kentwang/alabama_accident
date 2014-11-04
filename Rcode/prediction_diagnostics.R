load("./data/absError.RData")

#== check multiple boxplot
boxplot(absError)
summary(absError)


#== try median absolute error
score = data.frame(poisReg=medae(mu.cv.poisReg,Y),poisReg.offset=medae(mu.cv.poisReg.offset,Y),
                   negBino=medae(mu.cv.negBino,Y),negBino.offset=medae(mu.cv.negBino.offset,Y),
                   gbmPois=min(medae(mu.cv.gbmPois,Y)),gbmPois.offset=min(medae(mu.cv.gbmPois.offset,Y)),
                   tree=medae(mu.cv.tree,Y),tree.offset=medae(mu.cv.tree.offset,Y),
                   glmnet=min(medae(mu.cv.glmnet,Y)),glmnet.offset=min(medae(mu.cv.glmnet.offset,Y)))
print(round(score,4))

plot(1:length(score), score, "h", ylab = "MEDAE", xlab = "", 
     xaxt = "n", lwd = 2, main = "Comparison of model families on MEDAE")
points(1:length(score), score, pch = 19)
axis(1, at = 1:10, labels=FALSE)
text(x = 1:10, y = par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
     labels = names(score), srt = 30, adj = 1, xpd =TRUE)

