#############################################################################
### Mohamed Omar 
### 11/5/2019
### Goal: Creating a predictive model for bladder cancer progression using Logistic Regression
###         [Elastic Net Regression]
###############################################################################

## Clean Working environment
rm(list = ls())

## Set working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Logistic_Regression")

## Load necessary libraries
library(glmnet)
library(caret)
library(limma)
library(pROC)

############################################################

## Load data
load("./Objs/progressionDataGood.rda")

## Normalization between arrays
usedTrainMat <- normalizeBetweenArrays(mixTrainMat)
usedTestMat <- normalizeBetweenArrays(mixTestMat)

## The same for groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

##############################################################

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
names_train <- c(as.vector(rownames(usedTrainMat)))
colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))

### Fitting the model (Elastic Net Regression)
# Family ="binomial" means that this will be a logistic regression model
# Note: if we want linear regression, family will be ="gaussian"
# alpha = 0.5 means that this will be Elastic Net Regression model 
# Note: if we want ridge regression, alpha will be 0, AND if we want Lasso regression, alpha will be 1

set.seed(333)
fit <- cv.glmnet(Training, usedTrainGroup, type.measure = "class", alpha=0.5, family="binomial")

## Plot the fit model
png(filename = "./Figs/elasticNetRegression_model.png", width = 2000, height = 2000, res = 400)
plot(fit)
dev.off()

## Plot the fit model coefficints with gene labels
library(plotmo)
png("./Figs/coefficients_GeneLabels.png", width = 2000, height = 2000, res = 400)
plot_glmnet(fit$glmnet.fit)
dev.off()

####################################################################
## Access the names of predictors (genes) used to build the model
tmp_coeffs <- coef(fit, s=fit$lambda.1se)
Predictor_genes <- data.frame(name=tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1], coefficient = tmp_coeffs@x)
Predictor_genes_sorted <- Predictor_genes[order(Predictor_genes$coefficient, decreasing = TRUE),]
save(Predictor_genes_sorted, file = "./Objs/Predictor_genes_ENR.rda")

############################################################

## Plot the ROC curve and AUC in the training set

# Get predictions in the training set
train_pred_responses <- predict(fit, s=fit$lambda.1se, newx = Training, type="response")
# Plot the ROC curve
png("./Figs/ROC_Train.png", width = 2000, height = 2000, res = 300)
roc(usedTrainGroup, train_pred_responses[,1], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main="ELastic Net Regression ROC Train")
dev.off()

## Other method but gives different results !!
set.seed(333)
fit_auc <- cv.glmnet(Training, usedTrainGroup, type.measure = "auc", alpha=0.5, family="binomial")
AUC_train <- max(fit_auc$cvm)
AUC_train  ## 83% !! while it is 98.6% above 
#####################################################################

## Using the model to predict in the testing set 

# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
names_Test <- c(as.vector(rownames(usedTestMat)))
colnames(Testing) <- names_Test
names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

#### Plot the ROC and AUC curve in the Testing set
# Get predictions(responses) in the testing set
test_pred_Responses <- predict(fit, s=fit$lambda.1se, newx = Testing, type="response")
# Then plot the ROC curve
png("./Figs/ROC_Test.png", width = 2000, height = 2000, res = 300)
roc(usedTestGroup, test_pred_Responses[,1], plot = TRUE, print.auc=TRUE, levels = c("NoProgression", "Progression"), direction = ">", col="blue", lwd=2, grid=TRUE, main="ELastic Net Regression ROC Test")
dev.off()

## Class Predictions in the testing set
test_pred_classes <- predict(fit, s=fit$lambda.1se, newx = Testing, type="class")
table(test_pred_classes)
names_prediction <- rownames(test_pred_classes)
## Converting test_pred to factor instead of matrix
test_pred_classes <- as.vector(test_pred_classes)
names(test_pred_classes) <- names_prediction
test_pred_classes <- ordered(test_pred_classes, levels=c("Progression", "NoProgression"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test <- confusionMatrix(test_pred_classes, usedTestGroup, positive = "Progression")
Confusion_test