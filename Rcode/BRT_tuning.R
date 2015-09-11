gbmGrid = expand.grid(interaction.depth = c(1, 2, 3, 4),
	n.trees = seq(500,10000,by=200),
	shrinkage = c(0.01, 0.005, 0.001),
	n.minobsinnode = 20)

fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10)

gbmFit2 <- train(fmla.offset, data = a, method = "gbm",
	trControl = fitControl, tuneGrid = gbmGrid)

################################################
# use the old fashion in the paper
# gbm_6.offset_01 looks the best
################################################
gbm_3.offset_01 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=3,
                            shrinkage = .01)

gbm_3.offset_005 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=3,
                            shrinkage = .005)

gbm_3.offset_001 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=3,
                            shrinkage = .001)


gbm_4.offset_01 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=4,
                            shrinkage = .01)

gbm_4.offset_005 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=4,
                            shrinkage = .005)

gbm_4.offset_001 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=4,
                            shrinkage = .001)


gbm_5.offset_01 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=5,
                            shrinkage = .01)

gbm_5.offset_005 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=5,
                            shrinkage = .005)

gbm_5.offset_001 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=5,
                            shrinkage = .001)

gbm_6.offset_01 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=6,
                            shrinkage = .01)

gbm_6.offset_005 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=6,
                            shrinkage = .005)

gbm_6.offset_001 <- cv.gbmPois(fmla.offset, data = a, fold = fold,
                            tree.seq = tree.seq, interaction.depth=6,
                            shrinkage = .001)


plot(tree.seq,mlogL(gbm_3.offset_01,Y),type="l", ylim = c(1.5, 2.0), col=1,ylab="mlogL")
lines(tree.seq,mlogL(gbm_3.offset_005,Y), col=2)
lines(tree.seq,mlogL(gbm_3.offset_001,Y), col=3)
lines(tree.seq,mlogL(gbm_4.offset_01,Y), col=4)
lines(tree.seq,mlogL(gbm_4.offset_005,Y), col=5)
lines(tree.seq,mlogL(gbm_4.offset_001,Y), col=6)

title("GBM Comparison Using Mean Loglikelihood")
legend("topleft" ,
	legend = c("Dept 3, shrinkage 0.01", "Dept 3, shrinkage 0.005", "Dept 3, shrinkage 0.001",
			"Dept 4, shrinkage 0.01", "Dept 4, shrinkage 0.005", "Dept 4, shrinkage 0.001"),
	col = 1:6, lwd = 2, text.font = 3, cex = 0.8)


plot(tree.seq,mlogL(gbm_5.offset_01,Y),type="l", ylim = c(1.5, 2.0), col=1,ylab="mlogL")
lines(tree.seq,mlogL(gbm_5.offset_005,Y), col=2)
lines(tree.seq,mlogL(gbm_5.offset_001,Y), col=3)
lines(tree.seq,mlogL(gbm_6.offset_01,Y), col=4)
lines(tree.seq,mlogL(gbm_6.offset_005,Y), col=5)
lines(tree.seq,mlogL(gbm_6.offset_001,Y), col=6)

title("GBM Comparison Using Mean Loglikelihood")
legend("topleft" ,
	legend = c("Dept 5, shrinkage 0.01", "Dept 5, shrinkage 0.005", "Dept 5, shrinkage 0.001",
			"Dept 6, shrinkage 0.01", "Dept 6, shrinkage 0.005", "Dept 6, shrinkage 0.001"),
	col = 1:12, lwd = 2, text.font = 3, cex = 0.8)











