curdate = "08-03-2020"

library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(SuperLearner)

dir_out <- "./results/two_yr_mort/"
process_data_dir <- "./data/cache_data/"
date_of_SL_YD <- "08-03-2020"
SL_YD_PATH <- paste0(dir_out,"newlog-NES-", date_of_SL_YD, ".RData")
MED_OUTFILE <- paste0(dir_out,"mediation-newlog-NES", curdate, ".RData")

final_dt <- read.fst("./data/cache_data/final_hosp_mort_2000.fst",
                     as.data.table=T)

# EDIT 7/15/2020: Only include zip codes with > 2 individuals
# EDIT 7/15/2020: Only look at the first 5 years of entry ages (65-69)
final_dt <- final_dt[entry_age_break=="65to69"]
final_dt <- final_dt[N >= 5]
final_dt[, entry_age_break := NULL]


NORTHEAST=c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")  
SOUTH=c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
MIDWEST=c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
WEST=c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")

final_dt[,region := if_else(statecode %in% NORTHEAST, "NORTHEAST", 
                            ifelse(statecode %in% SOUTH, "SOUTH",
                                   ifelse(statecode %in% MIDWEST, "MIDWEST",
                                          ifelse(statecode %in% WEST, "WEST",
                                                 NA))))
         ]
final_dt <- na.omit(final_dt)


candidate.mediators <- c('1.1', '1.2', '1.3', '1.4', '1.5', '2.1', '2.2', 
                         '2.3', '2.4', '2.5', '2.6', '2.7', '2.8', '2.9', 
                         '2.10', '2.11', '2.12', '2.13', '2.14', '2.15', '2.16', 
                         '3.1', '3.2', '3.3', '3.4', '3.5', '3.6', '3.7', '3.8', 
                         '3.9', '3.10', '3.11', '4.1', '4.2', '4.3', '4.4', 
                         '5.1', '5.2', '5.3', '5.4', '5.5', '5.6', '5.7', '5.8',
                         '5.9', '5.10', '5.11', '5.12', '5.13', '5.14', '5.15', 
                         '6.1', '6.2', '6.3', '6.4', '6.5', '6.6', '6.7', '6.8',
                         '6.9', '7.1', '7.2', '7.3', '7.4', '7.5', '8.1', '8.2',
                         '8.3', '8.4', '8.5', '8.6', '8.7', '8.8', '8.9', '9.1',
                         '9.2', '9.3', '9.4', '9.5', '9.6', '9.7', '9.8', '9.9',
                         '9.10', '9.11', '9.12', '10.1', '10.2', '10.3', '11.1',
                         '11.2', '11.3', '11.4', '11.5', '11.6', '11.7', '12.1',
                         '12.2', '12.3', '12.4', '13.1', '13.2', '13.3', '13.4',
                         '13.5', '13.6', '13.7', '13.8', '13.9', '14.1', '14.2',
                         '14.3', '14.4', '14.5', '15.1', '15.2', '15.3', '15.4', 
                         '15.5', '15.6', '15.7', '16.1', '16.2', '16.3', '16.4', 
                         '16.5', '16.6', '16.7', '16.8', '16.9', '16.10', '16.11',
                         '16.12', '17.1', '17.2', "other_hosp")

treat.col <- "pm25_ensemble"
outcome.col <- "surv"
x.cols <- c('sex', 'race', 'dual', #'entry_age_break', 
            'statecode', 'mean_bmi', 
            'smoke_rate', 'hispanic', 'pct_blk', 'medhouseholdincome', 'medianhousevalue', 
            'poverty', 'education', 'popdensity', 'pct_owner_occ', 'summer_tmmx', 
            'winter_tmmx', 'summer_rmax', 'winter_rmax')
treat.pkey <- c("zip")
stratum.cols <- c("sex","race","dual"#,"entry_age_break"
)

msum=apply(final_dt[,..candidate.mediators], 2, sum)
zero_med_names <- names(msum[msum == 0])
candidate.mediators <- setdiff(candidate.mediators, zero_med_names)
final_dt <- final_dt[,-..zero_med_names]

# transform mediators into rates with YEAR_COUNT
final_M <- final_dt[, ..candidate.mediators]
# final_M <- final_dt[,
#                     lapply(.SD, function(x) x/final_dt$YEAR_COUNT),
#                     .SDcols = candidate.mediators
#                     ]
final_noM <- final_dt[,-..candidate.mediators]

# # rename mediators to begin with M
setnames(final_M, old=candidate.mediators, new=paste0("M", candidate.mediators))
candidate.mediators <- paste0("M", candidate.mediators)

# other preprocessing
# switch from survival to mortality
final_noM[, THREEYRSURV := N - THREEYRSURV]
final_noM[, surv := THREEYRSURV/N]
# final_noM[, pm25_ensemble := 1*(pm25_ensemble > 10)]
final_noM[, pm25_ensemble := 1*(pm25_ensemble > 13)]
final_noM[, c(stratum.cols) := lapply(.SD, as.factor),
          .SDcols = stratum.cols
          ]
final_noM[, zip := as.factor(zip)]
final_noM[, statecode := as.factor(statecode)]
final_noM[, region := as.factor(region)]

final_dt <- data.frame(final_noM, final_M) %>%
  filter(region %in% c("NORTHEAST","SOUTH"))

set.seed(234098)
# set.seed(3450981)
# set.seed(9813048)
all_zips <- unique(final_dt$zip)
selected_zip_len <- floor(length(all_zips)/2)
selected_zips <- all_zips[sample(length(all_zips), size=selected_zip_len)]

final_dt <- final_dt %>%
  filter(zip %in% selected_zips)


#### Now Remove Constant M from the subsample ####
msum=apply(final_dt[,candidate.mediators], 2, sum)
near_zero_med_names <- names(msum[msum <= 10])
candidate.mediators <- setdiff(candidate.mediators, near_zero_med_names)
final_dt <- final_dt %>%
  select(-one_of(near_zero_med_names))


surv_mean = mean(final_dt$surv)
omega=0.001

print(paste0("Sample size: ", nrow(final_dt)))

##### Poisson Analysis #####
library(gnm)
gnm_bin<-gnm(as.formula(
              paste0("THREEYRSURV ~ pm25_ensemble + ",
                    paste0(candidate.mediators, collapse=" + "), " + ",
                    paste0(x.cols, collapse=" + "), " + ",
                    "offset(log(N))"
                    )
              )
             ,
             data=final_dt,family=poisson(link="log"))
add_rate_coefs <- coef(gnm_bin); add_rate_coefs
Poisson_bin <- summary(gnm_bin); Poisson_bin
rm(gnm_bin)

Poisson_X <- model.matrix(
  as.formula(
    paste0(" ~ pm25_ensemble + ",
           paste0(x.cols, collapse=" + "), " + ",
           "offset(log(YEAR_COUNT))"
    )
  ), data=final_dt
)
Poisson_m <- sapply(
  candidate.mediators, function(med){
    mfit=glmnet::cv.glmnet(Poisson_X[,-1], pull(final_dt, med), family="poisson",
                           intercept=T)
    coef(mfit, s="lambda.min")[2]
#     mfit=glm(
#       as.formula(paste0("final_dt$", med, "~ Poisson_X - 1")),
#       family=poisson(link="log")
#     )
#     summary(mfit)$coefficients[2,]
  }
)
summary(Poisson_m)
# summary(Poisson_m[1,])

# # save(log_model, Poisson_bin_ne,file=paste0(dir_out,"log_model.RData"))

#### Robinson analysis
final_dt <- final_dt %>%
  mutate(surv = log((THREEYRSURV + omega*surv_mean)/(N + omega))) %>%
  mutate_at(vars(one_of(candidate.mediators)), 
            ~log((. + omega*mean(.))/(YEAR_COUNT + omega)))
# %>%
#   select(surv,
#          pm25_ensemble,
#          one_of(x.cols),
#          one_of(candidate.mediators),
#          zip
#          )

# linear confounding
gnm_log<-gnm(
              as.formula(
                paste0("surv ~ pm25_ensemble + ",
                      paste0(candidate.mediators, collapse=" + "), " + ",
                      paste0(x.cols, collapse=" + "))
              )
             ,
             #              eliminate= (as.factor(sex) + as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
             data=final_dt,family=gaussian())
log_lin <- summary(gnm_log)
rm(gnm_log)
print(log_lin)


#### Prepare for SuperLearner ####
rm(final_M, final_noM)

set.seed(234098)
folds <- caret::groupKFold(final_dt$zip, k=5) # returns train
folds <- lapply(folds, function(x) setdiff(1:nrow(final_dt), x))
lambdas <- 2^seq(-6, 6,length.out = 201L)

source("./code/Z_mediation_estimation_scripts.R")
source("./code/Z_caret_SL_def.R") # returns caret_learners, SL functions

start <- Sys.time()
tryCatch({
  
#   bin_lib <- c("SL.earth", "SL.glm", "SL.gam", "SL.mean", xgb_lrn$names, rf_bin_names)
#   cont_lib <- c("SL.earth", "SL.glm", "SL.gam", "SL.mean", xgb_lrn$names, rf_con_names)
  
  bin_lib <- c("SL.earth", "SL.glm", "SL.mean", rf_bin_names)
  cont_lib <- c("SL.earth", "SL.glm", "SL.mean", rf_con_names)
  
  set.seed(5460981)
  sl.mu.mx <- train_sl_mu(final_dt, candidate.mediators,
                          x.cols, treat.col, outcome.col,
                          #                        treat.pkey, stratum.cols,
                          bin_lib, cont_lib, folds, cores=5,
                          parallel_outfile="./cluster.log",
                          save_pred_only=TRUE,
                          onlyM=TRUE)
  save(sl.mu.mx,
       file=MED_OUTFILE)
  
  # load sl.mu for y and d
  load(SL_YD_PATH)
  
  # truncate
  sl.mu$mu.dx <- pmin(pmax(sl.mu$mu.dx, 0.05), 0.95)
  sl.mu <- c(sl.mu, sl.mu.mx)
  
  prd_result <- cv_sl_estimates(final_dt, mu.hats=sl.mu,
                                lambdas=lambdas,
                                folds=folds,
                                candidate.mediators=candidate.mediators, 
                                x.cols=x.cols, treat.col=treat.col, 
                                outcome.col=outcome.col,
                                weight.version="product")
  
  summary(prd_result$alpha_tildes)
  summary(prd_result$beta_tildes)
  summary(prd_result$alpha_hats * prd_result$beta_hats)
  prd_result$summ_Y
  prd_result$sel_M
  alpha_var <- sapply(
    prd_result$summ_allM[candidate.mediators %in% prd_result$sel_M], 
    function(x) vcov(x)[2,2])
  beta_cov <- vcov(prd_result$summ_Y)[-c(1:2),-c(1:2)]
  nde_var <- vcov(prd_result$summ_Y)[2,2]
  nie_var <- sum((prd_result$beta_hats^2)*alpha_var) + 
    sum((prd_result$alpha_hats^2)*diag(beta_var))
  nie_hat <- prd_result$NIE_hat
  nie_hat + c(-1,1)*1.96*sqrt(nie_var)
  nde_hat <- prd_result$NDE_hat
  nde_hat + c(-1,1)*1.96*sqrt(nde_var)
  
  mix_result <- cv_sl_estimates(final_dt, mu.hats=sl.mu,
                                lambdas=lambdas, 
                                folds=folds,
                                candidate.mediators=candidate.mediators, 
                                x.cols=x.cols, treat.col=treat.col, 
                                outcome.col=outcome.col,
                                weight.version="mixture")
  
  summary(mix_result$alpha_hats * mix_result$beta_hats)
  mix_result$sel_M
  alpha_var <- sapply(
    mix_result$summ_allM[candidate.mediators %in% mix_result$sel_M], 
    function(x) vcov(x)[2,2])
  beta_var <- vcov(mix_result$summ_Y)[-c(1:2),-c(1:2)]
  nde_var <- vcov(mix_result$summ_Y)[2,2]
  nie_var <- sum((mix_result$beta_hats^2)*alpha_var) + 
    sum((mix_result$alpha_hats^2)*diag(beta_var))
  nie_hat <- mix_result$NIE_hat
  nie_hat + c(-1,1)*1.96*sqrt(nie_var)
  nde_hat <- mix_result$NDE_hat
  nde_hat + c(-1,1)*1.96*sqrt(nde_var)
  
  foo=mix_result$alpha_hats * mix_result$beta_hats
  names(foo) <- mix_result$sel_M
  save(sl.mu, SL_YD_PATH, prd_result, mix_result,
       file=MED_OUTFILE)
}, 
error=function(e) e, 
finally={end <- Sys.time()}
)

print(paste0("Execution time: ", difftime(end, start)))

sl.mu$mu.dx.coef

sl.mu$mu.yx.coef
