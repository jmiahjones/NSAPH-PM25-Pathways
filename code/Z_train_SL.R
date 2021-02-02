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

train_sl_mu <- function(sel_df,
                        x.cols, treat.col, outcome.col,
                        bin_lib, cont_lib, folds, cores=NULL, 
                        parallel_outfile=NULL,
                        save_pred_only=TRUE,
                        d_method="method.NNloglik",
                        y_method="method.NNLS") {
  #require(foreach)
  V <- length(folds)
  SL.CV.control <- list(V=V, validRows=folds, shuffle=FALSE)
  saveAllFlag <- FALSE
  
  # selection data
  X <- model.matrix(as.formula("~."), 
                    data=sel_df %>% dplyr::select(x.cols))[,-1] %>%
    data.frame
  D <- dplyr::pull(sel_df, treat.col)
  
  # tests
  if(is.factor(dplyr::pull(sel_df, treat.col))){
    if(length(levels(dplyr::pull(sel_df, treat.col))) != 2){
      stop("Only accepting binary treatments!")
    }
    message("Assuming exposure is the second level.")
    active_lvl = levels(dplyr::pull(sel_df, treat.col))[2]
    D = 1*(dplyr::pull(sel_df, treat.col) == active_lvl)
    
  } else {
    lvls = sort(unique(dplyr::pull(sel_df, treat.col)))
    if(length(lvls) != 2) {
      stop("Only accepting binary treatments! Treatment does not have 2 unique values.")
    }
    if(!all.equal(lvls, c(0,1))){
    message("Assuming exposure is the largest value.")
    active_lvl = lvls[2]
    D = 1*(dplyr::pull(sel_df, treat.col) == active_lvl)
    }
  }
  
  
  parallelB = !is.null(cores) & is.numeric(cores)
  if(parallelB){
    require(parallel)
    if(!is.null(parallel_outfile)){
      if(file.exists(parallel_outfile)) 
        file.remove(parallel_outfile)
    }
    cl = parallel::makeCluster(cores, "FORK", outfile=parallel_outfile)
  } else {
    cl = "seq"
  }
  
  print("Estimating Mu.dx")
  
  mu.dx = try(
    CV.SuperLearner(Y=D, X=X, 
                    SL.library=bin_lib, family=binomial(),
                    cvControl=SL.CV.control,
                    saveAll=saveAllFlag,
                    parallel=cl,
                    method=d_method)
  )
  
  if(inherits(mu.dx, "try-error")){
    message("D Method failed: Trying backup.")
    mu.dx <- CV.SuperLearner(Y=D, X=X, 
                             SL.library=bin_lib, family=binomial(),
                             cvControl=SL.CV.control,
                             saveAll=saveAllFlag, parallel=cl,
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
  
  mu.yx = try(
              CV.SuperLearner(Y=dplyr::pull(sel_df, outcome.col), X=X,
                          SL.library=cont_lib, family=gaussian(),
                          cvControl=SL.CV.control,
                          saveAll=saveAllFlag, parallel=cl,
                          method=y_method)
  )
  
  if(inherits(mu.yx, "try-error")){
    message("Y Method failed: Trying backup.")
    mu.yx = try(
                CV.SuperLearner(Y=dplyr::pull(sel_df, outcome.col), X=X,
                          SL.library=cont_lib, family=gaussian(),
                          cvControl=SL.CV.control,
                          saveAll=saveAllFlag, parallel=cl,
                          method="method.CC_LS")
    )
  }
  if(inherits(mu.yx, "try-error")){
    stop("SuperLearner for Y|X Failed!")  
  }
  if(save_pred_only){
    mu.yx.pred = my_sl_predict(mu.yx)
    mu.yx.coef = mu.yx$coef
    rm(mu.yx)
  }
  
  if(parallelB){
    stopCluster(cl)
  }
  
  
  return(
    list(mu.dx=mu.dx.pred,
         mu.yx=mu.yx.pred,
         mu.dx.coef=mu.dx.coef,
         mu.yx.coef=mu.yx.coef)
  )
  
}


train_sl_mu_grouped_D <- function(sel_df,
                                  x.cols, treat.col, outcome.col,
                                  treat.pkey, stratum.cols,
                                  bin_lib, cont_lib, folds, cores=NULL, 
                                  parallel_outfile=NULL,
                                  save_pred_only=TRUE) {
  #require(foreach)
  V <- length(folds)
  SL.CV.control <- list(V=V, validRows=folds, shuffle=FALSE)
  saveAllFlag <- FALSE
  
  # selection data
  X <- sel_df %>% dplyr::select(x.cols)
  
  x.minus.strat.cols <- setdiff(x.cols, stratum.cols)
  
  trt_df <- dplyr::select(sel_df, all_of(c(treat.col, treat.pkey, x.minus.strat.cols))) %>%
    dplyr::distinct(.keep_all=T)
  X_trt <- dplyr::select(trt_df, all_of(c(x.minus.strat.cols)))
  D <- dplyr::pull(trt_df, treat.col)
  
  # tests
  if(is.factor(dplyr::pull(trt_df, treat.col))){
    if(length(levels(dplyr::pull(trt_df, treat.col))) != 2){
      stop("Only accepting binary treatments!")
    }
    message("Assuming exposure is the second level.")
    active_lvl = levels(dplyr::pull(trt_df, treat.col))[2]
    D = 1*(dplyr::pull(trt_df, treat.col) == active_lvl)
    
  } else {
    lvls = sort(unique(dplyr::pull(trt_df, treat.col)))
    if(length(lvls) != 2) {
      stop("Only accepting binary treatments! Treatment does not have 2 unique values.")
    }
    if(!all.equal(lvls, c(0,1))){
      message("Assuming exposure is the largest value.")
      active_lvl = lvls[2]
      D = 1*(dplyr::pull(trt_df, treat.col) == active_lvl)
    }
  }
  
  
  parallelB = !is.null(cores) & is.numeric(cores)
  if(parallelB){
    require(parallel)
    cl = parallel::makeCluster(cores, "FORK", outfile=parallel_outfile)
  } else {
    cl = "seq"
  }
  
  print("Estimating Mu.dx")
  
  mu.dx = try(
    CV.SuperLearner(Y=D, X=X_trt, 
                    SL.library=bin_lib, family=binomial(),
                    cvControl=SL.CV.control,
                    saveAll=saveAllFlag,
                    parallel=cl,
                    method="method.NNloglik")
  )
  
  if(inherits(mu.dx, "try-error")){
    mu.dx <- CV.SuperLearner(Y=D, X=X_trt, 
                             SL.library=bin_lib, family=binomial(),
                             cvControl=SL.CV.control,
                             saveAll=saveAllFlag, parallel=cl,
                             method="method.CC_LS")
  }
  if(inherits(mu.dx, "try-error")){
    stop("SuperLearner for D|X Failed!")  
  }
  if(save_pred_only){
    mu.dx = my_sl_predict(mu.dx)
  }
  
  joined_Dhat <- trt_df %>%
    mutate(pred=mu.dx) %>%
    dplyr::select(all_of(c(treat.pkey, "pred"))) %>%
    dplyr::right_join(sel_df, by=treat.pkey) %>%
    dplyr::pull(pred)
  rm(trt_df)
  
  if(length(joined_Dhat) != nrow(sel_df)){
    warning(paste0("Error joining treatment to outcome: ",
                   "treatment of length ", length(joined_Dhat),
                   " while outcome of length ", nrow(sel_df)))
  }
  
  print("Estimating Mu.yx")
  
  mu.yx = try(
    CV.SuperLearner(Y=dplyr::pull(sel_df, outcome.col), X=X,
                    SL.library=cont_lib, family=gaussian(),
                    cvControl=SL.CV.control,
                    saveAll=saveAllFlag, parallel=cl,
                    method="method.NNLS")
  )
  
  if(inherits(mu.yx, "try-error")){
    mu.yx = try(
      CV.SuperLearner(Y=dplyr::pull(sel_df, outcome.col), X=X,
                      SL.library=cont_lib, family=gaussian(),
                      cvControl=SL.CV.control,
                      saveAll=saveAllFlag, parallel=cl,
                      method="method.CC_LS")
    )
  }
  if(inherits(mu.yx, "try-error")){
    stop("SuperLearner for Y|X Failed!")  
  }
  if(save_pred_only){
    mu.yx = my_sl_predict(mu.yx)
  }
  
  if(parallelB){
    stopCluster(cl)
  }
  
  
  return(
    list(mu.dx=joined_Dhat,
         mu.yx=mu.yx)
  )
  
}
