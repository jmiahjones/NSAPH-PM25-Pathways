# go through the whole mediation analysis with gam; no selection
curdate = "09-08-2020"

library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(fst)
library(foreach)
library(doParallel)
# library(SuperLearner)

dir_out <- "./results/two_yr_mort/"
process_data_dir <- "./data/cache_data/"



dir_out <- "./results/two_yr_mort/"
process_data_dir <- "./data/cache_data/"
overlapdate <- "08-18-2020"
OVERLAP_PATH <- paste0(dir_out,"newlog-NES-", overlapdate, ".RData")
SL_YD_PATH <- paste0(dir_out,"gam-NES-", curdate, ".RData")
MED_OUTFILE <- paste0(dir_out,"mediation-gam-NES", curdate, ".RData")

final_dt <- read.fst("./data/cache_data/final_hosp_mort_2000.fst",
                     as.data.table=T)

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
final_noM[, pm25_ensemble := 1*(pm25_ensemble > 12)]
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


# get overlap variable
load(OVERLAP_PATH)
final_dt <- final_dt[overlap,]

#### Now Remove Constant M from the subsample ####
msum=apply(final_dt[,candidate.mediators], 2, sum)
near_zero_med_names <- names(msum[msum <= 10])
candidate.mediators <- setdiff(candidate.mediators, near_zero_med_names)
final_dt <- final_dt %>%
  select(-one_of(near_zero_med_names))


# surv_mean = mean(final_dt$surv)
# omega=0.001

print(paste0("Sample size: ", nrow(final_dt)))
# 
# ##### Binomial GAM Analysis #####
# library(gam)
# start <- Sys.time()
# gam_formula <- as.formula(paste("cbind(THREEYRSURV, I(N - THREEYRSURV))~ pm25_ensemble + ",
#                                 paste(candidate.mediators, collapse = " + "), " + ",
#                                 paste(paste("s(", 
#                                             x.cols[-c(1:4)], ",", 2, ")", sep = ""), 
#                                       collapse = " + "), 
#                                 "+", 
#                                 paste(x.cols[1:4], collapse = " + ")))
# 
# gam_y <- gam(gam_formula, data=final_dt, family=binomial())
# nde.hat <- coef(gam_y)[2]
# betas <- coef(gam_y)[(2 + (1:length(candidate.mediators)))]
# 
# cl <- parallel::makeCluster(5)
# registerDoParallel(cl)
# 
# alphas <- foreach(m=candidate.mediators, .combine=c,
#                   .export=c("final_dt", "x.cols"),
#                   .packages="gam") %dopar% {
#   gam_formula <- as.formula(paste("cbind(", m, ", YEAR_COUNT) ~ pm25_ensemble + ", 
#                                   paste(paste("s(", 
#                                               x.cols[-c(1:4)], ",", 2, ")", sep = ""), 
#                                         collapse = " + "), 
#                                   "+", 
#                                   paste(x.cols[1:4], collapse = " + ")))
#   
#   gam_m <- gam(gam_formula, data=final_dt, family=binomial())
#   coef(gam_m)[2]
# }
# stopImplicitCluster()
# nie.paths <- exp(alphas)*exp(betas)
# nie.hat <- sum(nie.paths)
# 
# end <- Sys.time()
# 
# print(paste0("Execution time: ", difftime(end, start)))
# difftime(end, start)
# 
# nde.hat
# nie.hat
# 
# gam_formula <- as.formula(paste("cbind(THREEYRSURV, I(N - THREEYRSURV))~ pm25_ensemble + ",
#                                 paste(paste("s(", 
#                                             x.cols[-c(1:4)], ",", 2, ")", sep = ""), 
#                                       collapse = " + "), 
#                                 "+", 
#                                 paste(x.cols[1:4], collapse = " + ")))
# 
# gam_y_nom <- gam(gam_formula, data=final_dt, family=binomial())
# coef(gam_y_nom)[2]


##### Poisson GAM Analysis #####
library(gam)
start <- Sys.time()

cl <- parallel::makeCluster(5)
registerDoParallel(cl)

alphas <- foreach(m=candidate.mediators, .combine=c,
                  .export=c("final_dt", "x.cols"),
                  .packages="gam") %dopar% {
                    gam_formula <- as.formula(paste(m, " ~ pm25_ensemble + ", 
                                                    paste(paste("s(", 
                                                                x.cols[-c(1:4)], ",", 2, ")", sep = ""), 
                                                          collapse = " + "), 
                                                    "+", 
                                                    paste(x.cols[1:4], collapse = " + "),
                                                    "+ offset(log(YEAR_COUNT))"))
                    
                    gam_m <- gam(gam_formula, data=final_dt, family=poisson())
                    coef(gam_m)[2]
                  }
stopImplicitCluster()

# m0s <- foreach(m.idx=seq_along(candidate.mediators), .combine=cbind) %do% {
#   m <- final_dt %>% dplyr::select(candidate.mediators[m.idx]) %>%
#     as.matrix
#   m0 <- m * exp(- (alphas[m.idx] * final_dt$pm25_ensemble))
#   return(m0)
# }

gam_formula <- as.formula(paste("THREEYRSURV~ pm25_ensemble + ",
                                paste(candidate.mediators, collapse=" + "), " + ",
                                paste(paste("s(", 
                                            x.cols[-c(1:4)], ",", 2, ")", sep = ""), 
                                      collapse = " + "), 
                                "+", 
                                paste(x.cols[1:4], collapse = " + "),
                                " + offset(log(N))"))

gam_y <- gam(gam_formula, data=final_dt, family=poisson())
betas <- coef(gam_y)[(2 + (1:length(candidate.mediators)))]
theta <- coef(gam_y)[2]
nde.hat <- exp(theta)
nie.paths <- exp(exp(alphas)*(betas))
nie.hat <- prod(nie.paths)

nde.hat
nie.hat
nie.hat*nde.hat

# mdiff <- foreach(m.idx=seq_along(candidate.mediators), .combine=cbind) %do% {
#   final_dt$pm25_ensemble * exp(alphas[m.idx])
# }
# colnames(mdiff) <- candidate.mediators
# 
# gam_formula <- as.formula(paste("THREEYRSURV~ pm25_ensemble + ",
#                                 " mdiff + ",
#                                 paste(paste("s(", 
#                                             x.cols[-c(1:4)], ",", 2, ")", sep = ""), 
#                                       collapse = " + "), 
#                                 "+", 
#                                 paste(x.cols[1:4], collapse = " + "),
#                                 " + offset(log(N))"))
# 
# gam_y <- gam(gam_formula, data=final_dt, family=poisson())
# nde.hat <- coef(gam_y)[2]
# betas <- coef(gam_y)[(2 + (1:length(candidate.mediators)))]
# nie.paths <- exp(exp(alphas)*(betas))
# nie.hat <- prod(nie.paths)



end <- Sys.time()

print(paste0("Execution time: ", difftime(end, start)))
difftime(end, start)

gam_formula <- as.formula(paste("THREEYRSURV~ pm25_ensemble + ",
                                paste(paste("s(", 
                                            x.cols[-c(1:4)], ",", 2, ")", sep = ""), 
                                      collapse = " + "), 
                                "+", 
                                paste(x.cols[1:4], collapse = " + "),
                                " + offset(log(N))"))

gam_y_nom <- gam(gam_formula, data=final_dt, family=poisson())
total.eff <- exp(coef(gam_y_nom)[2]); total.eff

