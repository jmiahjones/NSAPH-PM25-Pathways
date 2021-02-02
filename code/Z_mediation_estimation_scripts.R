# pseudo-mse-with-ml-estimation-scripts.R
# Functions that perform the estimation of NDE and NIE using semiparametric confounding control.
# Includes both the double-robust linear technique, as well as the svm method.

#################################################################
### Methods for computing cross-validation
#################################################################

##### Cross-validated SuperLearner #####

# shorter version
cv_sl_estimates <- function(sel_df, mu.hats,
                            weight.version, lambdas, folds,
                            candidate.mediators, x.cols, treat.col, outcome.col) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  require(dplyr)
  require(SuperLearner)
  
  
  p = length(candidate.mediators)
  V <- length(folds)
  
  mu.mxis <- mu.hats$mu.mxis
  mu.dx <- mu.hats$mu.dx
  mu.yx <- mu.hats$mu.yx
  
  # tests
  if(is.factor(dplyr::pull(sel_df, treat.col))){
    if(length(levels(dplyr::pull(sel_df, treat.col))) != 2){
      stop("Only accepting binary treatments!")
    }
    message("Assuming exposure is the second level.")
    active_lvl = levels(dplyr::pull(sel_df, treat.col))[2]
    sel_df[,treat.col] = 1*(dplyr::pull(sel_df, treat.col) == active_lvl)
    
  } else {
    lvls = sort(unique(dplyr::pull(sel_df, treat.col)))
    if(length(lvls) != 2) {
      stop("Only accepting binary treatments! Treatment does not have 2 unique values.")
    }
    if(!all.equal(lvls, c(0,1))){
    message("Assuming exposure is the largest value.")
    active_lvl = lvls[2]
    sel_df[,treat.col] = 1*(dplyr::pull(sel_df, treat.col) == active_lvl)
    }
  }
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  m_0 = sapply(1:p, function(i){
    # estimate the counfounding effect on mediatior i
    mi = candidate.mediators[i]
    Mi <- dplyr::pull(sel_df, mi)
    mu.mxi = mu.mxis[[i]]
    
    col <- Mi - my_sl_predict(mu.mxi)
    return(col)
  })
  stopifnot(ncol(m_0)==p)
  
  d.prob = as.numeric(my_sl_predict(mu.dx))
  dc = dplyr::pull(sel_df, treat.col) - d.prob
  
  mfit <- lm(m_0 ~ dc)
  alpha_tildes = coef(mfit)[2, ]
  summ_mfit <- summary(mfit)
  rm(mfit)
  
  y_0 = as.numeric(dplyr::pull(sel_df, outcome.col) - my_sl_predict(mu.yx))
  message(class(y_0))
  yfit = lm(y_0 ~ dc + m_0)
  beta_tildes = coef(yfit)[-(1:2)]
  
  
  # decide which weight to use
  switch (weight.version,
          "mixture" = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes))),
          "product" = weights <- 1/(abs(beta_tildes) * abs(alpha_tildes)),
          "adaptive"= weights <- 1/(abs(beta_tildes))
          # "both"    = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes)))
  )
  
  print("Running Adaptive Lasso")
  
  relaxo.loss <- function(x, y, intercept=FALSE, family, penalty.factor, lambda,
                          holdout.x, holdout.y){
    lasso.fit <- glmnet::glmnet(x, y, intercept=intercept,
                                family = family, 
                                penalty.factor = penalty.factor,
                                lambda = lambda)
    
    beta_0 <- rep(0, ncol(x)+1)
    lapply(lambda, function(lam){
      lasso.coef <- coef(lasso.fit, s=lam)
      is_selected <- lasso.coef[-1,1] != 0
      beta <- beta_0
      beta[c(T, is_selected)] <- coef(lm(y~x[,is_selected]))
      
      error <- holdout.y - beta[1] - as.numeric(holdout.x %*% matrix(beta[-1], ncol=1))
      loss <- error^2
      return(loss)
    }) %>% do.call(cbind, .)
  }
  
  design_mat <- cbind(dc, m_0)
  losses <- vector("list", V)
  for(fold.idx in 1:V){
    holdout.idx <- folds[[fold.idx]]
    ## perform variable selection
    # regress mediators w/ d,x removed on y w/ d,x removed
    holdout.mat <- design_mat[holdout.idx,]
    holdout.y <- y_0[holdout.idx]
    train.mat <- design_mat[-holdout.idx,]
    train.y <- y_0[-holdout.idx]
    losses[[fold.idx]] <- relaxo.loss(x=train.mat, y=train.y, intercept=FALSE, 
                                      family = "gaussian", 
                                      penalty.factor = c(0, weights), 
                                      lambda = lambdas, holdout.x=holdout.mat,
                                      holdout.y = holdout.y)
  }
  loss <- do.call(rbind, losses) %>% colMeans()
  
  min.loss.lambda <- lambdas[which.min(loss)]
  
  adapt.fit <- glmnet::glmnet(x=design_mat, y=y_0, intercept=F, 
                              family="gaussian", penalty.factor = c(0, weights),
                              lambda=min.loss.lambda*(2^(-1:3)))
  
  adapt.fit.coef <- coef(adapt.fit, s=min.loss.lambda)
  is_selected <- adapt.fit.coef[-(1:2),] != 0
  
  
  
  # print("Getting Estimated Models")
  
  if(!any(is_selected)){
    alpha_hats = NA
    beta_hats = NA
    NIE_hat = NA
    yfit = lm(y_0 ~ dc)
    gamma_hat = coef(yfit)[-1]
  } else {
    selected_idx <- which(is_selected)
    
    ## first stage after variable selection
    ## must use the estimation data set
    
    # center the new treatments
    m_0_selected = m_0[, selected_idx]
    # alpha_hats = coef(lm(m_0_selected ~ dc))[2, ]
    alpha_hats = alpha_tildes[selected_idx]
    
    yfit = lm(y_0 ~ dc + m_0_selected)
    gamma_hat <- coef(yfit)[2]
    beta_hats <- coef(yfit)[-(1:2)]
    
    NIE_hat <- sum(alpha_hats * beta_hats)
  }
  
  
  return(
    list(alpha_hats = alpha_hats, beta_hats = beta_hats,
         alpha_tildes=alpha_tildes, beta_tildes=beta_tildes,
         NIE_hat = NIE_hat, NDE_hat = gamma_hat,
         ATE_hat = NIE_hat + gamma_hat,
         sel_M = candidate.mediators[is_selected],
         losses=loss,
         lambda=min.loss.lambda,
         summ_Y = summary(yfit),
         summ_allM = summ_mfit)
  )
  
}

#####################################################
# Cross-validation: Compute estimates at each lambda
#####################################################
# shorter version
cv_sl_estimates_tuning <- function(sel_df, mu.hats,
                                   weight.version, lambdas, folds,
                                   candidate.mediators, x.cols, treat.col, outcome.col) {
  
  stopifnot(length(lambdas) > 1) # assumption. Need to change returns
  require(dplyr)
  require(SuperLearner)
  
  
  p = length(candidate.mediators)
  V <- length(folds)
  
  mu.mxis <- mu.hats$mu.mxis
  mu.dx <- mu.hats$mu.dx
  mu.yx <- mu.hats$mu.yx
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  m_0 = sapply(1:p, function(i){
    # estimate the counfounding effect on mediatior i
    mi = candidate.mediators[i]
    Mi <- dplyr::pull(sel_df, mi)
    mu.mxi = mu.mxis[[i]]
    
    col <- Mi - my_sl_predict(mu.mxi)
    return(col)
  })
  stopifnot(ncol(m_0)==p)
  
  d.prob = as.numeric(my_sl_predict(mu.dx))
  dc = sel_df[,treat.col] - d.prob
  
  alpha_tildes = coef(lm(m_0 ~ dc))[2, ]
  
  y_0 = sel_df[,outcome.col] - my_sl_predict(mu.yx)
  yfit = lm(y_0 ~ dc + m_0)
  beta_tildes = coef(yfit)[-(1:2)]
  
  
  # decide which weight to use
  switch (weight.version,
          "mixture" = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes))),
          "product" = weights <- 1/(abs(beta_tildes) * abs(alpha_tildes)),
          "adaptive"= weights <- 1/(abs(beta_tildes))
          # "both"    = weights <- 1/(abs(beta_tildes) * (1 + abs(alpha_tildes)))
  )
  
  print("Running Adaptive Lasso")
  
  relaxo.loss <- function(x, y, intercept=FALSE, family, penalty.factor, lambda,
                          holdout.x, holdout.y){
    lasso.fit <- glmnet::glmnet(x, y, intercept=intercept,
                                family = family, 
                                penalty.factor = penalty.factor,
                                lambda = lambda)
    
    beta_0 <- rep(0, ncol(x)+1)
    lapply(lambda, function(lam){
      lasso.coef <- coef(lasso.fit, s=lam)
      is_selected <- lasso.coef[-1,1] != 0
      beta <- beta_0
      relax.fit <- lm(y~x[,is_selected])
      epsilons <- residuals(relax.fit)
      beta[c(T, is_selected)] <- coef(relax.fit)
      error <- holdout.y - beta[1] - as.numeric(holdout.x %*% matrix(beta[-1], ncol=1))
      # loss <- error^2
      
      ## ASSUME: D is in the first column in X, so X[,-1] is the centered M
      # M <- X[,-1]
      loss <- sapply(1:p, function(j){
        abs(alpha_tildes[j] * sum(error*holdout.x[,j+1]))
      }) %>% sum
      return(loss)
    }) %>% do.call(cbind, .)
  }
  
  design_mat <- cbind(dc, m_0)
  losses <- vector("list", V)
  for(fold.idx in 1:V){
    holdout.idx <- folds[[fold.idx]]
    ## perform variable selection
    # regress mediators w/ d,x removed on y w/ d,x removed
    holdout.mat <- design_mat[holdout.idx,]
    holdout.y <- y_0[holdout.idx]
    train.mat <- design_mat[-holdout.idx,]
    train.y <- y_0[-holdout.idx]
    losses[[fold.idx]] <- relaxo.loss(x=train.mat, y=train.y, intercept=FALSE, 
                                      family = "gaussian", 
                                      penalty.factor = c(0, weights), 
                                      lambda = lambdas, holdout.x=holdout.mat,
                                      holdout.y = holdout.y)
  }
  loss <- do.call(rbind, losses) %>% colSums()
  
  min.loss.lambda <- lambdas[which.min(loss)]
  total_effect = coef(lm(y_0 ~ dc))[-1]
  
  results_per_lambda <- lapply(lambdas, function(lam){
    adapt.fit <- glmnet::glmnet(x=design_mat, y=y_0, intercept=F, 
                                family="gaussian", penalty.factor = c(0, weights),
                                lambda=lam*(2^(-1:2)))
    
    adapt.fit.coef <- coef(adapt.fit, s=lam)
    is_selected <- adapt.fit.coef[-(1:2),] != 0
    
    if(!any(is_selected)){
      alpha_hats = NA
      beta_hats = NA
      NIE_hat = NA
      gamma_hat = total_effect
    } else {
      selected_idx <- which(is_selected)
      
      ## first stage after variable selection
      ## must use the estimation data set
      
      # center the new treatments
      m_0_selected = m_0[, selected_idx]
      # alpha_hats = coef(lm(m_0_selected ~ dc))[2, ]
      alpha_hats = alpha_tildes[selected_idx]
      
      yfit = lm(y_0 ~ dc + m_0_selected)
      gamma_hat <- coef(yfit)[2]
      beta_hats <- coef(yfit)[-(1:2)]
      
      NIE_hat <- sum(alpha_hats * beta_hats)
    }
    
    
    return(
      list(alpha_hats = alpha_hats, beta_hats = beta_hats,
           NIE_hat = NIE_hat, NDE_hat = gamma_hat,
           sel_M = candidate.mediators[is_selected],
           lambda=min.loss.lambda)
    )
  })
  
  return(
    list(min.loss.lambda=min.loss.lambda,
         results=results_per_lambda,
         loss=loss
    )
  )
}


my_sl_predict <- function(obj, Y=NULL){
  require(SuperLearner)
  obj_class <- class(obj)
  if(!obj_class %in% c("SuperLearner", "CV.SuperLearner")){
    #stop("Internal error: Only superlearners allowed for obj.")
    # treat these as already-computed predictions for each observation
    out <- obj
  }
  
  if(obj_class == "SuperLearner"){
    if(all(coef(obj) == 0))
    {
      if(is.null(Y)){
        stop("All metalearner coefficients are zero!")
      } else {
        out <- rep(mean(Y), length(Y))
      }
    } else {
      out <- predict(obj, onlySL=TRUE)$pred
    }
  } else if(obj_class == "CV.SuperLearner") {
    out <- obj$SL.predict
  }
  return(out)
}

train_sl_mu <- function(sel_df, candidate.mediators, 
                        x.cols, treat.col, outcome.col,
                        bin_lib, cont_lib, folds, cores=NULL, 
                        parallel_outfile=NULL,
                        save_pred_only=TRUE,
                        onlyM=FALSE) {
  #require(foreach)
  p = length(candidate.mediators)
  V <- length(folds)
  SL.CV.control <- list(V=V, validRows=folds, shuffle=FALSE)
  saveAllFlag <- FALSE
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  # tests
  if(is.factor(dplyr::pull(sel_df, treat.col))){
    if(length(levels(dplyr::pull(sel_df, treat.col))) != 2){
      stop("Only accepting binary treatments!")
    }
    message("Assuming exposure is the second level.")
    active_lvl = levels(dplyr::pull(sel_df, treat.col))[2]
    sel_df[,treat.col] = 1*(dplyr::pull(sel_df, treat.col) == active_lvl)
    
  } else {
    lvls = sort(unique(dplyr::pull(sel_df, treat.col)))
    if(length(lvls) != 2) {
      stop("Only accepting binary treatments! Treatment does not have 2 unique values.")
    }
    if(!all.equal(lvls, c(0,1))){
    message("Assuming exposure is the largest value.")
    active_lvl = lvls[2]
    sel_df[,treat.col] = 1*(dplyr::pull(sel_df, treat.col) == active_lvl)
    }
  }
  
  print("Estimating Mu.mx")
  
  parallelB = !is.null(cores) & is.numeric(cores)
  if(parallelB){
    require(doParallel)
    require(foreach)
    cl = makeCluster(cores, "FORK", outfile=parallel_outfile)
    registerDoParallel(cl)
  }
  # mu.mxis = lapply(1:p, function(i){
  mu.mxis = foreach(i=1:p) %dopar% {
#   mu.mxis = foreach(i=1:p) %do% {
    mi = candidate.mediators[i]
    print(paste0(Sys.time(), ": Starting ", i, " of ", p))
    Mi <- dplyr::pull(sel_df, mi)
    
    mu.mxi <- try(
      CV.SuperLearner(Y=Mi, X=X, SL.library=cont_lib, family=gaussian(),
                      cvControl=SL.CV.control,
                      saveAll=saveAllFlag,
                      method="method.CC_LS")
    )
    
    if(inherits(mu.mxi, "try-error")){
      mu.mxi <- CV.SuperLearner(Y=Mi, X=X, SL.library=cont_lib, family=gaussian(),
                                cvControl=SL.CV.control,
                                saveAll=saveAllFlag,
                                method="method.NNLS")
      
      if(inherits(mu.mxi, "try-error")){
        warning(paste0("SuperLearner for M|X Failed: ", mi))
        return(NULL)
      }
    }
    if(save_pred_only){
      mu.mxi.pred = my_sl_predict(mu.mxi)
      mu.mxi.coef = mu.mxi$coef
      rm(mu.mxi)
    }
    return(
      list(mu.mxi=mu.mxi.pred,
           mu.mxi.coef=mu.mxi.coef)
    )
  }
  mu.mxi.preds = lapply(mu.mxis, function(x) x$mu.mxi)
  mu.mxi.coefs = lapply(mu.mxis, function(x) x$mu.mxi.coef)
  rm(mu.mxis)
  
  if(parallelB){
    stopCluster(cl)
  }
  
  if(onlyM){
    return(
      list(mu.mxis=mu.mxi.preds,
           mu.mxi.coefs=mu.mxi.coefs)
    )
  } else {
    print("Estimating Mu.dx")
    
    mu.dx = try(
      CV.SuperLearner(Y=dplyr::pull(sel_df, treat.col), X=X, 
                      SL.library=bin_lib, family=binomial(),
                      cvControl=SL.CV.control,
                      saveAll=saveAllFlag,
                      method="method.NNloglik")
    )
    
    if(inherits(mu.dx, "try-error")){
      mu.dx <- CV.SuperLearner(Y=dplyr::pull(sel_df, treat.col), X=X, 
                               SL.library=bin_lib, family=binomial(),
                               cvControl=SL.CV.control,
                               saveAll=saveAllFlag,
                               method="method.CC_LS")
    }
    if(inherits(mu.dx, "try-error")){
      stop("SuperLearner for D|X Failed!")  
    }
    if(save_pred_only){
      mu.dx.pred = my_sl_predict(mu.dx)
      mu.dx.coef = mu.dx$coef
      rm(mu.dx)
    }
    
    
    print("Estimating Mu.yx")
    
    mu.yx = CV.SuperLearner(Y=dplyr::pull(sel_df, outcome.col), X=X,
                            SL.library=cont_lib, family=gaussian(),
                            cvControl=SL.CV.control,
                            saveAll=saveAllFlag,
                            method="method.CC_LS")
    
    if(inherits(mu.yx, "try-error")){
      stop("SuperLearner for Y|X Failed!")  
    }
    if(save_pred_only){
      mu.yx.pred = my_sl_predict(mu.yx)
      mu.yx.coef = mu.yx$coef
      rm(mu.yx)
    }
    
    return(
      list(mu.mxis=mu.mxi.preds,
           mu.dx=mu.dx.pred,
           mu.yx=mu.yx.pred,
           mu.mxi.coefs=mu.mxi.coefs,
           mu.dx.coef=mu.dx.coef,
           mu.yx.coef=mu.yx.coef)
    )
  }
}
