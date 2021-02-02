library(fst)
# library(data.table)
library(parallel)
library(dplyr)

cache_file_dir <- "./data/cache_data/"
covariate_dir <- paste0(cache_file_dir, "mortality_covariates/")
aggregate_mortality_dir <- paste0(cache_file_dir, "aggregate_mortality/")

f <- list.files("./data/mortality/",
                pattern = "\\.fst",
                full.names = TRUE)

myvars <- c("year","zip","sex","race","age","dual","entry_age_break","statecode",
            "followup_year","followup_year_plus_one","dead","pm25_ensemble",
            "mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue",
            "poverty","education","popdensity", "pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")
# national_merged2016 <- rbindlist(lapply(f,
#                                         read_fst,
#                                         columns = myvars,
#                                         as.data.table=TRUE))


NORTHEAST=c("NY","MA","PA","RI","NH","ME","VT","CT","NJ")  
SOUTH=c("DC","VA","NC","WV","KY","SC","GA","FL","AL","TN","MS","AR","MD","DE","OK","TX","LA")
MIDWEST=c("OH","IN","MI","IA","MO","WI","MN","SD","ND","IL","KS","NE")
WEST=c("MT","CO","WY","ID","UT","NV","CA","OR","WA","AZ","NM")

for(file in f){
  this.year <- read_fst(file, columns=myvars)
  this.year$zip <- sprintf("%05d", this.year$zip)
  
  
  this.year$region=ifelse(this.year$statecode %in% NORTHEAST, "NORTHEAST", 
                          ifelse(this.year$statecode %in% SOUTH, "SOUTH",
                                 ifelse(this.year$statecode %in% MIDWEST, "MIDWEST",
                                        ifelse(this.year$statecode %in% WEST, "WEST",
                                               NA))))
  
  this.year <- na.omit(this.year)
  
  # Covariates
  covariates <- this.year %>%
    group_by(year, zip) %>%
    summarize_at(vars(pm25_ensemble:region), ~min(.)) %>%
    na.omit()
  
  covariates$year_fac <- as.factor(covariates$year)
  covariates$region <- as.factor(covariates$region)
  
  cov_year <- covariates$year[1]
  
  covariate_file <- paste0(covariate_dir, "mortality_covariates_", 
                           cov_year, ".fst")
  write.fst(covariates,path=covariate_file)
  rm(covariates)
  
  
  # Generate count data for each individual characteristics and follow-up year
  this.year$time_count<-this.year$followup_year_plus_one-this.year$followup_year
  
  this.grouped.year <- this.year %>%
    group_by(zip,year,sex,race,dual,entry_age_break,followup_year)
  
  dead_personyear <- this.grouped.year %>%
    summarize_at(vars(dead, time_count), ~mean(.))
  
  confounders <- this.grouped.year %>%
    summarize_at(vars(pm25_ensemble:region), ~min(.))
  
  rm(this.grouped.year)
  
  aggregate_data <- inner_join(dead_personyear,confounders,
                               by=c("zip","year","sex","race","dual",
                                    "entry_age_break","followup_year")) %>%
    na.omit
  
  agg_year <- aggregate_data$year[1]
  agg_file <- paste0(aggregate_mortality_dir, "aggregate_mortality_", 
                     agg_year, ".fst")
  write.fst(aggregate_data,path=agg_file)
  
  rm(dead_personyear, confounders, aggregate_data)
  
}


# # All types of statistical models implemented for the paper
# library(gnm)
# # library("parallel")
# # require(doParallel)
# # library(data.table)
# # library(fst)
# 
# 
# files <- list.files(aggregate_mortality_dir,
#                     pattern = "\\.fst",
#                     full.names = TRUE)
# 
# national_data <- rbindlist(lapply(files, read_fst))
# 
# 
# # # Cox Proportional Hazard
# # Cox_raw <- coxph(Surv(followup_year,followup_year_plus_one,dead)~pm25_ensemble + 
# #                    mean_bmi + smoke_rate + hispanic + pct_blk +
# #                    medhouseholdincome + medianhousevalue +
# #                    poverty + education + popdensity + pct_owner_occ +
# #                    summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
# #                    as.factor(year) + as.factor(region)
# #                  +strata(as.factor(entry_age_break))+strata(as.factor(sex))+strata(as.factor(race))+strata(as.factor(dual))
# #                  , data=national_merged2016,
# #                  ties = c("efron"),na.action = na.omit)
# # Cox <- summary(Cox_raw)
# # save(Cox ,file=paste0(dir_out,"Cox.RData"))
# 
# 
# # Cox-equvalent conditional Poisson Regression
# gnm_raw<-gnm(dead~  pm25_ensemble + 
#                mean_bmi + smoke_rate + hispanic + pct_blk +
#                medhouseholdincome + medianhousevalue +
#                poverty + education + popdensity + pct_owner_occ +
#                summer_tmmx + winter_tmmx + summer_rmax + winter_rmax +
#                as.factor(year) + as.factor(region)
#              +offset(log(time_count)),
#              eliminate= (as.factor(sex):as.factor(race):as.factor(dual):as.factor(entry_age_break):as.factor(followup_year)), 
#              data=national_data,family=poisson(link="log"))
# Poisson<-summary(gnm_raw)
# exp(10*Poisson$coefficients[1])
# # save(Poisson,file=paste0(dir_out,"Poisson.RData"))