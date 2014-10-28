#-------------------------------------------------------------------------------
# Variable selection under offset models including 1) Poisson, 2) Negative Bino-
# mial, 3) glmnet Poisson , 5) poisson trees,  5) gbm Poisson.
#
# Bascially replicate what's been done in model_selection.r
#-------------------------------------------------------------------------------

# Poisson regression
poisReg.offset <- step(glm(fmla.offset, data = a, family = "poisson"))

# Negative Binomial
negBino.offset <- step(glm.nb(fmla.offset, data = a))

# glmnet
glmnetPois.cv.offset <- cv.glmnet(model.matrix(fmla.offset, data = a)[,-1], 
                                  as.matrix(a[, "X5YrCrashCount"]), family = "poisson", offset = log(a$Traffic))
beta.offset <- as.matrix(coef(glmnetPois.cv.offset, s="lambda.1se"))
beta.offset.names <- rownames(beta.offset)
beta.offset <- as.numeric(beta.offset)
names(beta.offset) <- beta.offset.names

# Poisson regression trees
treePois.offset <- rpart(fmla.offset, data = a, method = "poisson")

# gbm boosting trees
gbmPois.offset <- gbm(fmla.offset, data = a, distribution = "poisson", 
                      n.trees = 10000, interaction.depth = 3)
summary_gbm <- summary(gbmPois.offset)

