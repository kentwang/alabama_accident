#-------------------------------------------------------------------------------
# Variable selection under offset models including 1) Poisson, 2) Negative Bino-
# mial, 3) glmnet Poisson , 5) poisson trees,  5) gbm Poisson.
#
# Bascially replicate what's been done in model_selection.r
# (1) This is using the traditional way
#-------------------------------------------------------------------------------

# Poisson regression
poisReg.offset <- step(glm(fmla.offset, data = a, family = "poisson"))

# Negative Binomial
negBino.offset <- step(glm.nb(fmla.offset, data = a))

# glmnet
glmnetPois.cv.offset <- cv.glmnet(model.matrix(fmla.offset, data = a)[,-1], 
                                  as.matrix(a[, "X5YrCrashCount"]), family = "poisson", offset = log(a$Traffic))
beta.glmnet.offset <- as.matrix(coef(glmnetPois.cv.offset, s="lambda.min"))
beta.glmnet.offset.names <- rownames(beta.glmnet.offset)
beta.glmnet.offset <- as.numeric(beta.offset)
names(beta.glmnet.offset) <- beta.glmnet.offset.names
beta.glmnet.offset <- beta.glmnet.offset[-1] # exclude intercept

# Poisson regression trees
treePois.offset <- rpart(fmla.offset, data = a, method = "poisson")

# gbm boosting trees
gbmPois.offset <- gbm(fmla.offset, data = a, distribution = "poisson", 
                      n.trees = 10000, interaction.depth = 3)
summary_gbm <- summary(gbmPois.offset)
beta.gbm.offset <- summary_gbm$rel.inf
names(beta.gbm.offset) <- summary_gbm$var


# clarify variable names. offset variable not included
# run model_selection fmla.string part first
allVariables <- allVariables[allVariables != 'log(Traffic)']

## - Compare variable importances, using t-statistics all that different
# Note: use lambda.min in glmnet
varImp_poisReg_offset <- varImpStandard(summary(poisReg.offset)$coefficients[-1, 3])
varImp_negBino_offset <- varImpStandard(summary(negBino.offset)$coefficients[-1, 3])
varImp_glmnetPois_offset <- varImpStandard(beta.glmnet.offset)
varImp_treePois_offset <- varImpStandard(treePois.offset$variable.importance)
varImp_gbmPois <- varImpStandard(beta.gbm.offset)

varImp <- as.vector(c(varImp_poisReg_offset, varImp_negBino_offset,
                      varImp_glmnetPois_offset, varImp_treePois_offset, 
                      varImp_gbmPois))

model.names <- c("Poisson Regression", "Negative Binomial",
                 "glmnet Poisson", "Regression Tree", "gbm Boosted Trees")

importances <- as.data.frame(cbind(rep(allVariables, length(model.names)),
                                   rep(model.names, each = length(allVariables))))
importances$importance <- varImp
names(importances) <- c("var", "model", "importance")

ggplot(importances, aes(factor(var), importance, fill = model)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(text = element_text(size = 16)) + 
  coord_flip() + 
  scale_x_discrete(name="Intersection Factors") +
  scale_fill_discrete(breaks=model.names)

## - try use original coefficients instead of t-statistics
varImp_poisReg_offset <- varImpStandard(summary(poisReg.offset)$coefficients[-1, 3])
varImp_negBino_offset <- varImpStandard(summary(negBino.offset)$coefficients[-1, 3])
varImp_glmnetPois_offset <- varImpStandard(beta.glmnet.offset)

varImp <- as.vector(c(varImp_poisReg_offset, varImp_negBino_offset,
                      varImp_glmnetPois_offset))

model.names <- c("Poisson Regression", "Negative Binomial",
                 "glmnet Poisson")

importances <- as.data.frame(cbind(rep(allVariables, length(model.names)),
                                   rep(model.names, each = length(allVariables))))
importances$importance <- varImp
names(importances) <- c("var", "model", "importance")

ggplot(importances, aes(factor(var), importance, fill = model)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme(text = element_text(size = 16)) + 
  coord_flip() + 
  scale_x_discrete(name="Intersection Factors") +
  ggtitle("Variable importance in original scale")


################################################################################
# (2) Variable selection using the unified framework
# Note: the AIC calculation should be reversed
################################################################################
varImp.pois <- varImp.loo(fmla.offset, a, "poisson")
