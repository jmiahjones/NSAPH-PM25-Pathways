
library(NSAPHutils)
# set_threads()
library(icd)
library(fst)
library(data.table)
library(lubridate)
library(dtplyr)
library(dplyr)

dod_save_file <- "./data/cache_data/dod_2000_benefs_thru_2010.fst"

myvars <- c("qid","entry_year","zip","sex","race","dual","entry_age_break",
            "statecode","pm25_ensemble",
            "mean_bmi","smoke_rate","hispanic","pct_blk","medhouseholdincome","medianhousevalue",
            "poverty","education","popdensity", "pct_owner_occ","summer_tmmx","winter_tmmx","summer_rmax","winter_rmax")

# denom <- read_data("./data/mortality/", years = 2000, 
#                    columns = c("qid","zip","entry_year"))
denom <- read_data("./data/mortality/", years = 2000, 
                   columns = myvars)
# write down the entry age and race labels
denom[,entry_age_break:=factor(entry_age_break,
                               labels=c("65to69","70to74",
                                        "75to79","80to84",
                                        "85to89","90to94",
                                        "95over"))
      ]
# race labels, using D_MEDICARE_RACE coding (modified using RTI algorithm) 
# Source as of 2020: https://www.resdac.org/cms-data/variables/dmedicarerace
denom[,race:=factor(race, levels=c(0:6),
                    labels=c("Unknown",
                             "White",
                             "Black",
                             "Other",
                             "ASPAI",
                             "Hispanic",
                             "AmIndianAlaskNative"))]


# coarsen race categories
denom[,race:=factor(race, levels=c("White",
                                   "Black",
                                   "Unknown",
                                   "Other",
                                   "ASPAI",
                                   "Hispanic",
                                   "AmIndianAlaskNative"),
                    labels=c("White", "Black", "OtherUnknown", 
                             "OtherUnknown","OtherUnknown",
                             "OtherUnknown","OtherUnknown"))]


# only keep QIDs for those in the current entry year
denom <- denom[entry_year == 2000,][,entry_year := NULL]
year_2k_qids <- denom$qid
# denom <- denom %>% select(-entry_year)
# zip_lookup <- denom[,.(qid, zip)]
stratum_dt <- denom[,qid:entry_age_break]
confounder_dt <- denom[,-c("qid")]

##### Obtain Date of Death Info #####
# adm_file <- paste0("./data/admissions/admissions_by_year/admissions_", 
#                    2000, ".fst")
# dod_files <- paste0("./data/mortality/confounder_exposure_merged_nodups_health_", 
#                     2000L:2010L, ".fst")
# 
# out <- NULL
# for(dod_file in dod_files){
#   full_dod <- fst(dod_file)
#   
#   chunk_size <- 1000000L
#   start <- 1
#   end <- min(chunk_size, nrow(full_dod))
#   continue <- T
#   while (continue) {
#     temp <- read.fst(dod_file, from=start, to=end, 
# #                      columns=c("QID", "BENE_DOD"),
#                      columns=c("qid", "bene_dod"),
#                      as.data.table=T)
#     temp[,QID:=qid]
#     temp[,qid:=NULL]
#     temp[,BENE_DOD:=bene_dod]
#     temp[,bene_dod:=NULL]
#     temp <- temp[QID %in% year_2k_qids]
#     
#     if(nrow(temp) > 0){
#       temp[,BENE_DOD := dmy(BENE_DOD)]
#       temp <- temp[,lapply(.SD, max, na.rm=TRUE), by=.(QID)]
#       out <- rbind(out, temp)
#     }
#     
#     if (end >= nrow(full_dod)) {
#       continue <- F
#     }
#     
#     start <- end + 1
#     end <- min(end + chunk_size, nrow(full_dod))
#   }
#   # collapse the info from this year in with that of previous years
#   out <- out[,lapply(.SD, max, na.rm=TRUE), by=.(QID)]
# }
# rm(temp)
# 
# out <- out[QID %in% year_2k_qids]
# write.fst(out, path=dod_save_file)

# after first run, simply load the already-saved out file
out <- read.fst(dod_save_file, as.data.table=T)

## TODO: Change variable names
out[, DOD_YEAR := year(BENE_DOD)] # extract year, avoids weird Date NAs that aren't <NA>
out[, ALIVE := is.na(DOD_YEAR)]
out[, THREEYRSURV := 1*(BENE_DOD >= dmy("01JAN2002"))]
out[ALIVE==T, THREEYRSURV := 1]
out[, YEAR_COUNT := as.numeric(as.duration(
  difftime(dmy("01JAN2002"), dmy("01JAN2000"))
  ))/31536000
  ]
out[THREEYRSURV == 0, 
    YEAR_COUNT := as.numeric(as.duration(
      difftime(BENE_DOD, dmy("01JAN2000"))
      ))/31536000
    ]

summary(out$YEAR_COUNT)
table(out$THREEYRSURV)
out[THREEYRSURV == 1]
prop.table(table(out$THREEYRSURV))

hosp_count_file <- "./data/cache_data/hosps_by_qid/medicare_hosp_counts_entry_2000.fst"
foo=fst(hosp_count_file)
dim(foo)
hosp_counts <- read.fst(hosp_count_file, as.data.table=T)
unique(hosp_counts$AYEAR)
hosp_counts <- hosp_counts[AYEAR >= 2000 & AYEAR < 2002] # CHANGE THIS FOR DIFFERENT SURV
agg_hosp_counts <- hosp_counts[,
                               lapply(.SD, sum),
                               by = .(QID),
                               .SDcols=colnames(hosp_counts)[-(1:2)]
                               ]
# add in those without apparent hospitalizations
temp_qid <- setdiff(year_2k_qids, agg_hosp_counts$QID)

individual_count_df = out %>% 
  as_tibble %>%
  select(QID, THREEYRSURV, YEAR_COUNT) %>%
  left_join(as_tibble(agg_hosp_counts), by="QID") %>% 
  mutate_at(vars(-QID, -THREEYRSURV), 
            ~ifelse(is.na(.), 0, .)) %>% 
  inner_join(as_tibble(stratum_dt), by=c("QID"="qid")) %>%
  rename(other_hosp = ` `)

stratum_sum <- as.data.table(individual_count_df)[,-c("QID")][,
  lapply(.SD, sum),
  by = .(zip,sex,race,dual,entry_age_break)
  ]
#   lazy_dt(immutable=F) %>%
#   group_by(zip,sex,race,dual,entry_age_break) %>%
#   summarize_at(vars(THREEYRSURV:other_hosp), ~sum(.)) %>%
#   as_tibble

# stratum_denominator <- individual_count_df %>%
#   group_by(zip,sex,race,dual,entry_age_break) %>%
#   count
stratum_denominator <- as.data.table(individual_count_df)[,-c("QID")][,
  .N,
  by = .(zip,sex,race,dual,entry_age_break)
]

# confounder_df <- confounder_df %>%
#   lazy_dt %>%
#   group_by(zip,sex,race,dual,entry_age_break) %>%
#   summarize_all(~min(.)) %>%
#   as_tibble()
confounder_df <- confounder_dt[,
  lapply(.SD, min),
  by = .(zip,sex,race,dual,entry_age_break)
  ] %>% na.omit

final_dt <- merge(
    stratum_sum, 
    stratum_denominator,
    by=c("zip","sex","race","dual","entry_age_break"),
  ) %>%
  merge(
    confounder_df, 
    by=c("zip","sex","race","dual","entry_age_break")    
  )
final_dt[, zip := sprintf("%05d", zip)]
final_dt[, dual := as.integer(dual)]

final_dt <- final_dt[YEAR_COUNT > 0]

write_fst(final_dt, "./data/cache_data/final_hosp_mort_2000.fst")

# final_dt[,race:=factor(race, levels=c("White",
#                                       "Black",
#                                       "Unknown",
#                                       "Other",
#                                       "ASPAI",
#                                       "Hispanic",
#                                       "AmIndianAlaskNative"),
#                        labels=c("White", "Black", "OtherUnknown", 
#                                 "OtherUnknown","OtherUnknown","OtherUnknown","OtherUnknown"))]
# 
# coarse_dt <- final_dt[entry_age_break == "65to69" || entry_age_break == "70to74",
#                       lapply(.SD, sum), by = .(zip,sex,race,dual), 
#                       .SDcols=c("THREEYRSURV","N")]
# # coarse_dt <- coarse_dt[N > 5,] # around 7,000 observations
# nr <- nrow(coarse_dt)
# coarse_dt <- merge(coarse_dt,
#                    final_dt[,#entry_age_break == "65to69", 
#                             !c("entry_age_break","N","THREEYRSURV")
#                             ][,
#                               lapply(.SD, min),
#                               by=c("zip","sex","race","dual"),
#                             ],
#                    by=c("zip","sex","race","dual"))
# stopifnot(nr == nrow(coarse_dt))
# 
# write_fst(coarse_dt, "./data/cache_data/coarse_hosp_mort_2000.fst")
