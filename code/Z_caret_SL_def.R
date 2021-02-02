# my.SL.caret <- function (Y, X, newX, family, obsWeights, id, method = "rf", tuneLength = 5, 
#                          metric = ifelse(family$family == 
#                                            "gaussian", "RMSE", "Accuracy"), ...) 
# {
# #   .SL.require("caret")
#   
#   innerV <- 5
#   innerFold <- caret::groupKFold(id, k=innerV) # returns train
#   seeds <- sample.int(1e7, size=10*innerV*tuneLength)
#   trControl <- caret::trainControl(method = "cv", number = innerV,
#                                    index=innerFold,
#                                    verboseIter = FALSE,
#                                    trim=TRUE,
#                                    returnData=FALSE,
#                                    returnResamp="none",
#                                    seeds=seeds)
#   
#   if (family$family == "gaussian") {
#     fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
#                               metric = metric, method = method, tuneLength = tuneLength, 
#                               trControl = trControl)
#     pred <- predict(fit.train, newdata = newX, type = "raw")
#   }
#   if (family$family == "binomial") {
#     Y.f <- as.factor(Y)
#     levels(Y.f) <- c("A0", "A1")
#     fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
#                               metric = metric, method = method, tuneLength = tuneLength, 
#                               trControl = trControl)
#     pred <- predict(fit.train, newdata = newX, type = "prob")[, 
#                                                               2]
#   }
#   fit <- list(object = fit.train)
#   out <- list(pred = pred, fit = fit)
#   class(out$fit) <- c("SL.caret")
#   return(out)
# }

# caret_learners=SuperLearner::create.Learner("my.SL.caret", 
#                                             tune=list(method=c("rf", "nnet", "xgboost")),
#                                             detailed_names=T, name_prefix="SL.tuned")

rf_lrn <- create.Learner("SL.randomForest", 
                         tune=list(nodesize=c(2,5,10,15,20)),
                         detailed_names=T)
rf_bin_names <- rf_lrn$names[2]
rf_con_names <- rf_lrn$names[3:4]

xgb_lrn <- create.Learner("SL.xgboost", tune=list(ntrees=c(100)),
                          params=list(max_depth=2),
                          detailed_names=T)
