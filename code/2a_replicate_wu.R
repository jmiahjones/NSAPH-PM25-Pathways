library(data.table)
library(dplyr)
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


aggregate_data <- dplyr::left_join(aggregate_data,covariates,
                                   by=c("zip", "year"))

rm(covariates)

# Cox-equvalent conditional Poisson Regression
library(gnm)
gnm_raw<-gnm(dead~  pm25_ensemble + 
               mean_bmi + smoke_rate + hispanic + pct_blk +
               medhouseholdincome + medianhousevalue +
               poverty + education + popdensity + pct_owner_occ +
               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
               as.factor(year) + as.factor(region)
             +offset(log(time_count)),
             eliminate= (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
             data=aggregate_data,family=poisson(link="log"))
Poisson<-summary(gnm_raw)
rm(gnm_raw)
exp(10*Poisson$coefficients[1])
save(Poisson,file=paste0(dir_out,"Poisson.RData"))

gnm_bin<-gnm(dead ~ I(pm25_ensemble > 10) + 
               mean_bmi + smoke_rate + hispanic + pct_blk +
               medhouseholdincome + medianhousevalue +
               poverty + education + popdensity + pct_owner_occ +
               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
               as.factor(year) + as.factor(region)
             +offset(log(time_count)),
             eliminate= (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
             data=aggregate_data,family=poisson(link="log"))
exp(coef(gnm_bin)[1])
Poisson_bin_txt <- summary(gnm_bin)
rm(gnm_bin)
save(Poisson, Poisson_bin_txt,file=paste0(dir_out,"Poisson.RData"))

gnm_bin_region<-gnm(dead ~ I(pm25_ensemble > 10) + 
               mean_bmi + smoke_rate + hispanic + pct_blk +
               medhouseholdincome + medianhousevalue +
               poverty + education + popdensity + pct_owner_occ +
               summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
               as.factor(year) + as.factor(region) +
               I(pm25_ensemble > 10)*as.factor(region)
             +offset(log(time_count)),
             eliminate= (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
             data=aggregate_data,family=poisson(link="log"),
             start=c(Poisson_bin_txt$coefficients[,1], rep(NA, 3)))
Poisson_bin_region <- summary(gnm_bin_region)
exp(c(Poisson_bin_region$coefficients[1,1], 
      tail(Poisson_bin_region$coefficients[,1], n=3)))
rm(gnm_bin_region)
save(Poisson_bin_region,file=paste0(dir_out,"Poisson-region.RData"))