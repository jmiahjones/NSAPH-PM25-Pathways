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
#          region, time_count, dead,
         sex, race, dual,
         entry_age_break, followup_year, zip) %>%
  as.data.frame()


rm(covariates)

set.seed(234098)
# set.seed(3450981)
# set.seed(9813048)
all_zips <- unique(aggregate_data$zip)
selected_zips <- all_zips[sample(length(all_zips), size=400L)]

aggregate_data <- aggregate_data %>%
  filter(zip %in% selected_zips)

treat.col <- "pm25_ensemble"
outcome.col <- "log_mort"
x.cols <- setdiff(colnames(aggregate_data), c(treat.col, outcome.col, "zip"))

set.seed(234098)
folds <- caret::groupKFold(aggregate_data$zip, k=5) # returns train
folds <- lapply(folds, function(x) setdiff(1:nrow(aggregate_data), x))

source("./code/Z_train_SL.R")
source("./code/Z_caret_SL_def.R") # returns caret_learners, SL functions

bin_lib <- c("SL.earth", "SL.glm", "SL.gam", "SL.mean", rf_bin_names, xgb_lrn$names)
cont_lib <- c("SL.earth", "SL.glm", "SL.gam", "SL.mean", rf_con_names, xgb_lrn$names)

start <- Sys.time()

set.seed(5460981)
sl.mu <- train_sl_mu(aggregate_data,
                     x.cols, treat.col, outcome.col,
                     bin_lib, cont_lib, folds, cores=NULL, 
                     parallel_outfile="./cluster.log",
                     save_pred_only=TRUE)

dc <- aggregate_data$pm25_ensemble - sl.mu$mu.dx
yc <- aggregate_data$log_mort - sl.mu$mu.yx

lm.coef <- coef(lm(yc ~ dc))[2]

end <- Sys.time()

save(sl.mu,file=paste0(dir_out,"sl_full_subsample.RData"))

print(paste0("Execution time: ", difftime(end, start)))
print(lm.coef)
