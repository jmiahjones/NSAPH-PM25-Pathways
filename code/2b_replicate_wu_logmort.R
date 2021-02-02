library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(SuperLearner)
# f <- list.files("./data/cache_data/aggregate_mortality",
#                 pattern = "\\.fst",
#                 full.names = TRUE)
# national_data <- rbindlist(lapply(f,
#                                   read_fst,
#                                   as.data.table=TRUE))

dir_out <- "./results/wu_replication/"
process_data_dir <- "./data/cache_data/"
load(paste0(process_data_dir, "aggregate_data.RData"))
load(paste0(process_data_dir, "covariates.RData"))

covariates <- as.data.table(covariates)

aggregate_data <- aggregate_data %>%
  as.data.table %>%
  lazy_dt(immutable=FALSE) %>%
  filter(region == "NORTHEAST") %>%
  filter(followup_year <= 5, year <= 2005) %>%
  mutate(log_mort = log(dead/time_count + 0.01),
         pm25_ensemble = 1*(pm25_ensemble > 10),
         region = as.factor(region),
         sex = as.factor(sex),
         race = as.factor(race),
         dual = as.factor(dual),
         entry_age_break = as.factor(entry_age_break),
         followup_year = as.factor(followup_year)) %>%
  left_join(covariates, by=c("zip", "year")) %>%
  select(log_mort,
         pm25_ensemble,
         mean_bmi, smoke_rate, hispanic, pct_blk,
         medhouseholdincome, medianhousevalue,
         poverty, education, popdensity, pct_owner_occ,
         summer_tmmx, winter_tmmx, summer_rmax, winter_rmax,
         year, 
         region, time_count, dead,
         sex, race, dual,
         entry_age_break, followup_year, zip) %>%
  as.data.frame()

rm(covariates)

set.seed(234098)
# set.seed(3450981)
# set.seed(9813048)
all_zips <- unique(aggregate_data$zip)
selected_zips <- all_zips[sample(length(all_zips), size=200L)]

aggregate_data <- aggregate_data %>%
  filter(zip %in% selected_zips)

print(paste0("Sample size: ", nrow(aggregate_data)))
# library(gnm)
# start <- Sys.time()
# 
# 
# gnm_bin<-gnm(dead ~ pm25_ensemble 
# #              + mean_bmi + smoke_rate + hispanic + pct_blk +
# #                medhouseholdincome + medianhousevalue +
# #                poverty + education + popdensity + pct_owner_occ +
# #                summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
# #                as.factor(year)
#              +offset(log(time_count)),
#              eliminate= (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
#              data=aggregate_data,family=poisson(link="log"))
# pois_coefs <- coef(gnm_bin)
# # Poisson_bin_ne <- summary(gnm_bin)
# rm(gnm_bin)
# 
# 
# gnm_log<-gnm(log_mort ~ pm25_ensemble 
# #              + mean_bmi + smoke_rate + hispanic + pct_blk +
# #                medhouseholdincome + medianhousevalue +
# #                poverty + education + popdensity + pct_owner_occ +
# #                summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
# #                as.factor(year)
#              ,
#              eliminate= ((sex):(race):(dual):(entry_age_break):(followup_year)), 
#              data=aggregate_data,family=gaussian())
# log_coefs <- coef(gnm_log)
# # log_model <- summary(gnm_log)
# rm(gnm_log)
# end <- Sys.time()
# # save(log_model, Poisson_bin_ne,file=paste0(dir_out,"log_model.RData"))

treat.col <- "pm25_ensemble"
outcome.col <- "log_mort"
x.cols <- setdiff(colnames(aggregate_data), c(treat.col, outcome.col))

set.seed(234098)
folds <- caret::groupKFold(aggregate_data$zip, k=5) # returns train
folds <- lapply(folds, function(x) setdiff(1:nrow(aggregate_data), x))

source("./code/Z_train_SL.R")
source("./code/Z_caret_SL_def.R") # returns caret_learners, SL functions

start <- Sys.time()
tryCatch({
    
    bin_lib <- c("SL.earth", "SL.glm", "SL.gam", "SL.mean", xgb_lrn$names, rf_bin_names)
    cont_lib <- c("SL.earth", "SL.glm", "SL.gam", "SL.mean", xgb_lrn$names, rf_con_names)
    
    set.seed(5460981)
    sl.mu <- train_sl_mu(aggregate_data,
                         x.cols, treat.col, outcome.col,
                         bin_lib, cont_lib, folds, cores=5, 
                         parallel_outfile="./cluster.log",
                         save_pred_only=TRUE)
    
    dc <- aggregate_data$pm25_ensemble - sl.mu$mu.dx
    yc <- aggregate_data$log_mort - sl.mu$mu.yx
    # mu.dx <- sl.mu$mu.dx
    # mu.yx <- sl.mu$mu.yx
    
    lm.fit <- lm(yc ~ dc)
    lm.coef <- coef(lm.fit)[2]
    print(lm.coef)
  }, 
  error=function(e) e, 
  finally={end <- Sys.time()}
)


save(sl.mu, yc, dc,file=paste0(dir_out,"sl_full_subsample.RData"))

print(paste0("Execution time: ", difftime(end, start)))
