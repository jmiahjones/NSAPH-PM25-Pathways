curdate = "07-03-2020"

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
  #   filter(region == "NORTHEAST") %>%
  filter(followup_year <= 2, year <= 2002) %>%
  mutate(mort = dead/time_count,
         pm25_ensemble = 1*(pm25_ensemble > 10),
         region = as.factor(region),
         sex = as.factor(sex),
         race = as.factor(race),
         dual = as.factor(dual),
         entry_age_break = as.factor(entry_age_break),
         followup_year = as.factor(followup_year)) %>%
  left_join(covariates, by=c("zip", "year")) %>%
  as.data.frame()
# 
# rm(covariates)
# 
# set.seed(234098)
# # set.seed(3450981)
# # set.seed(9813048)
# all_zips <- unique(aggregate_data$zip)
# selected_zips <- all_zips[sample(length(all_zips), size=200L)]
# 
# aggregate_data <- aggregate_data %>%
#   filter(zip %in% selected_zips)

mort_mean = mean(aggregate_data$mort)
omega=0.001

#### log-transformed linear analysis
aggregate_data <- aggregate_data %>%
  mutate(mort = log((dead + omega*mort_mean)/(time_count + omega))) %>%
  select(mort,
         pm25_ensemble,
         mean_bmi, smoke_rate, hispanic, pct_blk,
         medhouseholdincome, medianhousevalue,
         poverty, education, popdensity, pct_owner_occ,
         summer_tmmx, winter_tmmx, summer_rmax, winter_rmax,
         year, 
         #          region, time_count, dead,
         sex, race, dual,
         entry_age_break, followup_year, zip)

# linear confounding
gnm_log<-gnm(mort ~ pm25_ensemble
             + mean_bmi + smoke_rate + hispanic + pct_blk +
               medhouseholdincome + medianhousevalue +
               poverty + education + popdensity + pct_owner_occ +
               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
               as.factor(year)
             + (as.factor(sex) + as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year))
             ,
             #              eliminate= (as.factor(sex) + as.factor(race)+as.factor(dual)+as.factor(entry_age_break)+as.factor(followup_year)), 
             data=aggregate_data,family=gaussian())
log_lin <- summary(gnm_log)
rm(gnm_log)
print(log_lin)

save(log_lin, file=paste0(dir_out,"newlog_linear_full_",curdate,".Rdata"))