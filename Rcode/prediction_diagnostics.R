fmla <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
                    LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
                    LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                    MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                    NumLanes + Offset + OffsetDist + OneWay + PaveType + 
                    PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
                    RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                    log(Traffic) + TurnProhib + Lat + Long") 

fmla.offset <- as.formula("X5YrCrashCount ~ AreaType + IntCat + IntTCType + LegRtType + 
                    LegSpeed + LegTCType + LegType + LegWidth + Lighting + 
                    LTLanes + LTLnLength + LTOffset + LTWidth + MedType + 
                    MedWidth + MergeLanes + NextPIDist + NumberLegs + 
                    NumLanes + Offset + OffsetDist + OneWay + PaveType + 
                    PedCross + RTChannel + RTLanes + RTLnLength + RTMoveCtrl + 
                    RTWidth + Rumble + SightLt + SightRt + SkewAngle + Terrain + 
                    offset(log(Traffic)) + TurnProhib + Lat + Long") 


Y = a$X5YrCrashCount  # response
pred = matrix(, 200, 2)

for (i in 1:200) {
  cat('iter ', i, '\n')
  fold = cvfolds(nrow(a),k=5)  # get cv partition info fold 20 and 5
  
  train = which(fold != 5)
  test = which(fold == 5)
  
  fit1 <- glm.my(fmla, data = data[train, ], family = poisson)
  fit2 <- glm.my(fmla.offset, data = data[train, ], family = poisson)
  
  pred1 <- predict(fit1, newdata = data[test, ], type = "response")
  pred2 <- predict(fit2, newdata = data[test, ], type = "response")
  
  pred[i, ] = cbind(mae(pred1, Y[test]), mae(pred2, Y[test]))
}

pred = cbind(pred1[pred1 < 10], pred2[pred2 < 10])
boxplot(pred1, pred2)
apply(pred, 2, median)

#== check multiple boxplot
score = data.frame(poisReg=abs(mu.cv.poisReg - Y), poisReg.offset=abs(mu.cv.poisReg.offset - Y), 
              negBino=abs(mu.cv.negBino - Y), negBino.offset=abs(mu.cv.negBino.offset - Y), 
              gbmPois=abs(mu.cv.gbmPois - Y), gbmPois.offset=abs(mu.cv.gbmPois.offset - Y), 
              tree=abs(mu.cv.tree - Y), tree.offset=abs(mu.cv.tree.offset - Y), 
              glmnet=abs(mu.cv.glmnet[, 1] - Y), glmnet.offset=abs(mu.cv.glmnet.offset[, 1] - Y))
# boxplot(score, ylim = c(0, 2))
summary(score)

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

#== check the distribution of the absolute error
absError = abs(mu.cv.poisReg - Y)
hist(absError, 1000)
summary(absError)
