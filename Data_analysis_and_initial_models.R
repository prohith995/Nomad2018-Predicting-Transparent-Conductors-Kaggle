setwd("C:/Users/user/Desktop/OR/ISLR/Assignment 3")
library(fpp2)
library(glmnet)
library(DAAG)
train = read.csv("train.csv", header = TRUE)
test  = read.csv("test.csv",  header = TRUE)

View(test)
dim(train)
dim(test)

#Logistic
logistic_model  = glm(SeriousDlqin2yrs ~ . ,  data = train, family = binomial)
logistic_model$converged

predictions = predict(logistic_model, newdata = test, type = "response" )

my_logistic_submission = data.frame(id = 1:101503, 	Probability = predictions)
write.csv(my_logistic_submission, "Logistic regression1.csv", row.names = F)

#NaiveBayes
library(e1071)
nb_model       = naiveBayes(SeriousDlqin2yrs ~ ., data = train)
nb_y       = predict(nb_model,       newdata = test, type = "raw")[,2]

my_nb_submission = data.frame(id = 1:101503, 	Probability = nb_y)
write.csv(my_nb_submission, "Naive Bayes1.csv", row.names = F)

#Lasso and Ridge
y_train = train$SeriousDlqin2yrs
set.seed(10^6)
lasso_model  = cv.glmnet( x = as.matrix(train[,-1]), y = y_train, alpha = 1, type.measure = "auc", 
                          family = "binomial")
ridge_model  = cv.glmnet( x = as.matrix(train[,-1]), y = y_train, alpha = 0, type.measure = "auc", 
                          family = "binomial")

lasso_y    = predict(lasso_model, newx = as.matrix(test), s = lasso_model$lambda.min, type = "response")
ridge_y    = predict(ridge_model, newx = as.matrix(test), s = ridge_model$lambda.min, type = "response")

my_lasso_submission = data.frame(id = 1:101503, 	Probability = lasso_y)
write.csv(my_lasso_submission, "Lasso Regression1.csv", row.names = F)

my_ridge_submission = data.frame(id = 1:101503, 	Probability = ridge_y)
write.csv(my_ridge_submission, "Ridge Regression1.csv", row.names = F)


#Tree model
library(rpart)

tree.model = rpart(SeriousDlqin2yrs ~ ., data = train)
tree_predictions = predict(tree.model, newdata = test)

my_tree_submission = data.frame(id = 1:101503, 	Probability = tree_predictions)
write.csv(my_tree_submission, "Tree model1.csv", row.names = F)


#random forests
library(randomForest)
fit <- randomForest(as.factor(SeriousDlqin2yrs) ~ ., data=train,importance=TRUE,ntree=2000)
x.ct.pred <- predict(fit, newdata=test)

my_rf_submission = data.frame(id = 1:101503, 	Probability = x.ct.pred)
write.csv(my_rf_submission, "RandomForest.csv", row.names = F)

#Bagging
library(ipred)
x.ip <- bagging( SeriousDlqin2yrs~ ., data=train)
x.ip.prob <- predict(x.ip, type="prob", newdata=test)

my_bs_submission = data.frame(id = 1:101503, 	Probability = x.ip.prob)
write.csv(my_bs_submission, "Bootstrapping.csv", row.names = F)
