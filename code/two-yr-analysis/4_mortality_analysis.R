curdate = "08-18-2020"

library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(fst)
library(SuperLearner)

dir_out <- "./results/two_yr_mort/"
process_data_dir <- "./data/cache_data/"
SL_YD_PATH <- paste0(dir_out,"newlog-NES-", curdate, ".RData")

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

# msum=apply(final_dt[,..candidate.mediators], 2, sum)
# zero_med_names <- names(msum[msum == 0])
# candidate.mediators <- setdiff(candidate.mediators, zero_med_names)
# final_dt <- final_dt[,-..zero_med_names]
# 
# # transform mediators into rates with YEAR_COUNT
# final_M <- final_dt[,
#                     lapply(.SD, function(x) x/final_dt$YEAR_COUNT),
#                     .SDcols = candidate.mediators
#                     ]
final_noM <- final_dt[,-..candidate.mediators]

# # rename mediators to begin with M
# setnames(final_M, old=candidate.mediators, new=paste0("M", candidate.mediators))
# candidate.mediators <- paste0("M", candidate.mediators)

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

final_dt <- data.frame(final_noM) %>%
  filter(region %in% c("NORTHEAST","SOUTH"))



set.seed(234098)
# set.seed(3450981)
# set.seed(9813048)
all_zips <- unique(final_dt$zip)
selected_zip_len <- floor(length(all_zips)/2)
selected_zips <- all_zips[sample(length(all_zips), size=selected_zip_len)]

final_dt <- final_dt %>%
  filter(zip %in% selected_zips)

surv_mean = mean(final_dt$surv)
omega=0.001

print(paste0("Sample size: ", nrow(final_dt)))

##### Overlap Population #####

set.seed(234097)
folds <- caret::groupKFold(final_dt$zip, k=3) # returns train
folds <- lapply(folds, function(x) setdiff(1:nrow(final_dt), x))

source("./code/Z_train_SL.R")
source("./code/Z_caret_SL_def.R") # returns caret_learners, SL functions


#   bin_lib <- c("SL.earth", "SL.glm", "SL.gam", "SL.mean", xgb_lrn$names, rf_bin_names)
#   cont_lib <- c("SL.earth", "SL.glm", "SL.gam", "SL.mean", xgb_lrn$names, rf_con_names)

bin_lib <- c("SL.earth", "SL.glm", "SL.mean", rf_bin_names)
cont_lib <- c("SL.earth", "SL.glm", "SL.mean", rf_con_names)

cl <- parallel::makeCluster(3, "FORK", outfile="./cluster.log")
overlap.d.sl =
  CV.SuperLearner(Y=final_dt$pm25_ensemble, 
                  X=model.matrix(~., dplyr::select(final_dt, x.cols))[,-1] %>% 
                    data.frame,
                  SL.library=bin_lib, family=binomial(),
                  cvControl=list(V=length(folds), 
                                 validRows=folds,
                                 shuffle=FALSE),
                  saveAll=F, parallel=cl,
                  method="method.NNloglik")
parallel::stopCluster(cl)
message("Propensities trained. Finding overlap.")

overlap.coef <- coef(overlap.d.sl)
overlap.mu.dx = my_sl_predict(overlap.d.sl)
rm(overlap.d.sl)
overlap <- overlap.mu.dx >= .05 & overlap.mu.dx <= .95
final_dt <- final_dt[overlap,]

print(paste0("Overlap Sample size: ", nrow(final_dt)))

##### Poisson Analysis #####
library(gnm)
gnm_bin<-gnm(
  as.formula(
    paste0("THREEYRSURV ~ pm25_ensemble + ",
           paste0(x.cols, collapse=" + "), " + ",
           "offset(log(N))"
    )
  )
#   , eliminate = (sex:race:dual)
  ,
  data=final_dt,family=poisson(link="log"))
add_rate_coefs <- coef(gnm_bin); add_rate_coefs
Poisson_bin <- summary(gnm_bin); Poisson_bin
rm(gnm_bin)

##### Binom Analysis #####
library(gnm)
gnm_binom<-gnm(
  as.formula(
    paste0("cbind(THREEYRSURV, I(N - THREEYRSURV)) ~ pm25_ensemble + ",
           paste0(x.cols, collapse=" + ")
    )
  )
#   , eliminate = (sex:race:dual)
  ,
  data=final_dt,family=binomial())
# add_rate_coefs <- coef(gnm_bin); add_rate_coefs
Binomial_summ <- summary(gnm_binom); Binomial_summ
rm(gnm_binom)

# # save(log_model, Poisson_bin_ne,file=paste0(dir_out,"log_model.RData"))

#### Robinson analysis
final_dt <- final_dt %>%
  mutate(surv = ((THREEYRSURV + omega*surv_mean)/(N + omega))) %>%
#   rate model
  mutate(surv=log(surv))
#   logistic model
#   mutate(surv=log(surv)/(log(1-surv)))
# %>%
#   select(surv,
#          pm25_ensemble,
#          one_of(x.cols),
#          zip
#   )

# linear confounding
gnm_log<-gnm(
  as.formula(
    paste0("surv ~ pm25_ensemble + ",
           paste0(x.cols, collapse=" + "))
  )
  ,
  #              eliminate= (as.factor(sex) + as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
  data=final_dt,family=gaussian())
log_lin <- summary(gnm_log)
rm(gnm_log)
print(log_lin)

set.seed(5460981)
newfolds <- caret::groupKFold(final_dt$zip, k=3) # returns train
newfolds <- lapply(newfolds, function(x) setdiff(1:nrow(final_dt), x))

start <- Sys.time()
tryCatch({
  
  sl.mu <- train_sl_mu(final_dt,
                       x.cols, treat.col, outcome.col,
                       #                        treat.pkey, stratum.cols,
                       bin_lib, cont_lib, newfolds, cores=3, 
                       parallel_outfile="./cluster.log",
                       save_pred_only=TRUE,
                       d_method="method.NNloglik",
                       y_method="method.NNLS")
  
  trunc.mu.dx <- pmin(pmax(sl.mu$mu.dx, 0.05), 0.95)
  dc <- final_dt$pm25_ensemble - trunc.mu.dx
  yc <- final_dt$surv - sl.mu$mu.yx
  # mu.dx <- sl.mu$mu.dx
  # mu.yx <- sl.mu$mu.yx
  
  lm.fit <- lm(yc ~ dc)
  lm.coef <- coef(lm.fit)[2]
  print(lm.coef)
  print(summary(lm.fit))
  
  save(add_rate_coefs, log_lin, sl.mu, yc, dc, overlap,
       file=SL_YD_PATH)
}, 
error=function(e) e, 
finally={end <- Sys.time()}
)

print(paste0("Execution time: ", difftime(end, start)))

print("Comparison:")
Poisson_bin$coefficients["pm25_ensemble",c(1,2,4)] %>% round(3)
Binomial_summ$coefficients["pm25_ensemble",c(1,2,4)] %>% round(3)
log_lin$coefficients["pm25_ensemble",c(1,2,4)] %>% round(3)
summary(lm.fit)$coefficients["dc",c(1,2,4)] %>% round(3)

ate_sl <- coef(lm.fit)[2]
set.seed(581071)
boot.idxs <- replicate(1000, sample(length(yc), replace=T), simplify=F)
boot_samp <- sapply(1:1000, function(b) 
  coef(lm(yc[boot.idxs[[b]]] ~ dc[boot.idxs[[b]]]))[2]
)
plot(density(boot_samp - ate_sl, kernel="gaussian"))
quantile(boot_samp, probs = c(0.025, 0.975))
sd(boot_samp - ate_sl)


#### using truncated propensities
poisson_propreg <- gnm(
  as.formula(
    paste0("THREEYRSURV ~ pm25_ensemble + ",
           "trunc.mu.dx + ",
           "offset(log(N))"
    )
  ),
  data=final_dt,family=poisson(link="log"))
coef(poisson_propreg)

binom_propreg <- gnm(
  as.formula(
    paste0("cbind(THREEYRSURV, I(N-THREEYRSURV)) ~ pm25_ensemble + ",
           "trunc.mu.dx"
    )
  ),
  data=final_dt,family=binomial())
coef(binom_propreg)

log_propreg <- gnm(
  as.formula(
    paste0("surv ~ pm25_ensemble + ",
           "trunc.mu.dx")
  )
  ,
  #              eliminate= (as.factor(sex) + as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
  data=final_dt,family=gaussian())
coef(log_propreg)

ipwY1 <- with(final_dt,
     mean(pm25_ensemble*mean(pm25_ensemble)/trunc.mu.dx*THREEYRSURV/N)
)
ipwY0 <- with(final_dt,
              mean((1-pm25_ensemble)*mean(1-pm25_ensemble)/(1-trunc.mu.dx)*THREEYRSURV/N)
)
log(ipwY1/ipwY0)

ipwY1 <- with(final_dt,
              mean(pm25_ensemble/trunc.mu.dx*THREEYRSURV/N)
)
ipwY0 <- with(final_dt,
              mean((1-pm25_ensemble)/(1-trunc.mu.dx)*THREEYRSURV/N)
)
log(ipwY1/ipwY0)

ipwY1 <- with(final_dt,
              sum(pm25_ensemble*1/trunc.mu.dx*THREEYRSURV/N)/
                sum(pm25_ensemble/trunc.mu.dx)
)
ipwY0 <- with(final_dt,
              sum((1-pm25_ensemble)*1/(1-trunc.mu.dx)*THREEYRSURV/N)/
                sum((1-pm25_ensemble)/(1-trunc.mu.dx))
)
log(ipwY1/ipwY0)

# #### using the overlap population
# # overlap <- sl.mu$mu.dx >= .05 & sl.mu$mu.dx <= .95
# ipwY1 <- with(final_dt[overlap,],
#               mean(pm25_ensemble*mean(pm25_ensemble)/sl.mu$mu.dx[overlap]*THREEYRSURV/N)
# )
# ipwY0 <- with(final_dt[overlap,],
#               mean((1-pm25_ensemble)*mean(1-pm25_ensemble)/(1-sl.mu$mu.dx[overlap])*THREEYRSURV/N)
# )
# log(ipwY1/ipwY0)
# 
# ipwY1 <- with(final_dt[overlap,],
#               mean(pm25_ensemble*1/sl.mu$mu.dx[overlap]*THREEYRSURV/N)
# )
# ipwY0 <- with(final_dt[overlap,],
#               mean((1-pm25_ensemble)*1/(1-sl.mu$mu.dx[overlap])*THREEYRSURV/N)
# )
# log(ipwY1/ipwY0)
# 
# 
# ipwY1 <- with(final_dt[overlap,],
#               sum(pm25_ensemble*1/sl.mu$mu.dx[overlap]*THREEYRSURV/N)/
#                 sum(pm25_ensemble/sl.mu$mu.dx[overlap])
# )
# ipwY0 <- with(final_dt[overlap,],
#               sum((1-pm25_ensemble)*1/(1-sl.mu$mu.dx[overlap])*THREEYRSURV/N)/
#                 sum((1-pm25_ensemble)/(1-sl.mu$mu.dx[overlap]))
# )
# log(ipwY1/ipwY0)
# 
# poisson_propreg <- gnm(
#   as.formula(
#     paste0("THREEYRSURV ~ pm25_ensemble + ",
#            "sl.mu$mu.dx[overlap] + ",
#            "offset(log(N))"
#     )
#   ),
#   data=final_dt[overlap,],family=poisson(link="log"))
# coef(poisson_propreg)
# 
# binom_propreg <- gnm(
#   as.formula(
#     paste0("cbind(THREEYRSURV, I(N-THREEYRSURV)) ~ pm25_ensemble + ",
#            "sl.mu$mu.dx[overlap]"
#     )
#   ),
#   data=final_dt[overlap,],family=binomial())
# coef(binom_propreg)
# 
# log_propreg <- gnm(
#   as.formula(
#     paste0("surv ~ pm25_ensemble + ",
#            "sl.mu$mu.dx[overlap]")
#   )
#   ,
#   #              eliminate= (as.factor(sex) + as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
#   data=final_dt[overlap,],family=gaussian())
# coef(log_propreg)
# 
# # the robinson model
# lm(yc[overlap] ~ dc[overlap])

# the gam models on whole / overlap population
library(gam)
gam_formula <- as.formula(paste("cbind(THREEYRSURV, I(N - THREEYRSURV))~ pm25_ensemble + ", 
                                paste(paste("s(", 
                                  x.cols[-c(1:4)], ",", 2, ")", sep = ""), 
                                  collapse = " + "), 
                                "+", 
                                paste(x.cols[1:4], collapse = " + ")))

gam_bin <- gam(gam_formula, data=final_dt, family=binomial())
coef(gam_bin)

# gam_bin_overlap <- gam(gam_formula, data=final_dt[overlap,], family=binomial())
# coef(gam_bin_overlap)


#### now, re-train superlearner on the overlap population ####
# file.remove("./cluster.log")
# set.seed(489719761)
# start <- Sys.time()
# newfolds <- caret::groupKFold(final_dt$zip[overlap], k=5) # returns train
# newfolds <- lapply(newfolds, function(x) setdiff(1:sum(overlap), x))
# cl <- parallel::makeCluster(3, "FORK", outfile="./cluster.log")
# mu.yx =
#   CV.SuperLearner(Y=final_dt$surv[overlap], 
#                   X=model.matrix(~., dplyr::select(final_dt[overlap,], x.cols))[,-1] %>% 
#                     data.frame,
#                   SL.library=c(cont_lib, xgb_lrn$names), family=gaussian(),
#                   cvControl=list(V=length(newfolds), 
#                                  validRows=newfolds,
#                                  shuffle=FALSE),
#                   saveAll=F, parallel=cl,
#                   method="method.NNLS")
# stopCluster(cl)
# end <- Sys.time()
# difftime(end, start)
# 
# mu.yx.pred = my_sl_predict(mu.yx)
# mu.yx.coef = mu.yx$coef
# rm(mu.yx)
# 
# yc.overlap <- final_dt$surv[overlap] - mu.yx.pred
# dc.overlap <- dc[overlap]
# lm(yc.overlap ~ dc.overlap)
