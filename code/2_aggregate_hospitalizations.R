## Create list of hospitalizations by QID

library(NSAPHutils)
set_threads()
library(icd)
library(fst)
library(data.table)
library(lubridate)

sharded_assign_multi_ccs <- function(hosps, icd_code, lvl=2) {
  ## Assumes as.icd has already been run
  if (all(names(hosps) != c("visit_num", "DIAG1"))) {
    stop("Hospitalization data does not seem properply prepared")
  } 
  out <- NULL
  start <- 1
  end <- min(1000000, nrow(hosps))
  while (TRUE) {
    temp <- hosps[start:end]
    if (icd_code == 9) {
      codes <- icd9_comorbid_ccs(temp, return_df = T, single=F, lvl=lvl)
    } else if (icd_code == 10) {
      codes <- icd10_comorbid_ccs(temp, return_df = T, single=F, lvl=lvl)
#       codes["0"] <- F #Code 0 removed in icd10, marking F, to allow combination with icd9
    } else {
      stop("icd_code must be either 9 or 10")
    }
#     names(codes)[2:285] <- paste0("ccs_", names(codes)[2:285])
    out <- rbind(out, codes)
    if (end == nrow(hosps)) {
      break()
    }
    start <- end + 1
    end <- min(end + 1000000, nrow(hosps))
  }
  return(out)
}

read_and_filter_hosps <- function(load_years, included_qids, columns){
  filenames <- paste0("./data/admissions/admissions_by_year/admissions_", 
                      load_years, ".fst")
  rbindlist(lapply(filenames, function(x) 
    read.fst(x, columns=columns, as.data.table=T)[QID %in% included_qids])
  )
}

start <- Sys.time()

# years <- 2000:2016
years <- 2008:2016
for (this_entry_year in years) {
  denom <- read_data("./data/mortality/", years = this_entry_year, 
                     columns = c("qid","zip","entry_year"))
  
  # only keep QIDs for those in the current entry year
  denom <- denom[entry_year == this_entry_year,]
  
  hosp_start_year <- max(2000, this_entry_year-2)
  hosp_end_year <- min(2016, this_entry_year+7)
  
  hosps <- read_and_filter_hosps(load_years=hosp_start_year:hosp_end_year, 
                                 included_qids=denom$qid, 
                                 columns = c("QID", "DIAG1", "ADATE", "DDATE"))
  hosps[, ADATE := dmy(ADATE)]
  hosps[, AYEAR := year(ADATE)]
  hosps[, ADATE := NULL]
  hosps[, DDATE := dmy(DDATE)]
  
#   hosps <- hosps[AYEAR == this_entry_year]
#   hosps <- merge(hosps, denom, by.x = "QID", by.y = "qid")
  hosps[,visit_num := 1:.N]
  
  ## DISCHARGES 10/1 AND AFTER USE ICD10, Other wise ICD 9
  codes <- NULL
  hosps_icd9 <- hosps[DDATE < "2015-10-01"]
  if(nrow(hosps_icd9) > 0){
    hosps_icd9[,DIAG1 := as.icd9(DIAG1)]
    codes <- rbind(codes,
                   sharded_assign_multi_ccs(hosps_icd9[,.(visit_num, DIAG1)], icd_code = 9)
    )
  }
  rm(hosps_icd9)
  hosps_icd10 <- hosps[DDATE >= "2015-10-01"]
  if(nrow(hosps_icd10) > 0){
    hosps_icd10[, DIAG1 := as.icd10(DIAG1)]
    codes <- rbind(codes,
                   sharded_assign_multi_ccs(hosps_icd10[,.(visit_num, DIAG1)], icd_code = 10)
    )
  }
  rm(hosps_icd10)
  if(is.null(codes) || nrow(codes) < 1){
    message(paste0("No hospitalizations found for year ", this_entry_year))
    warning(paste0("Skipping ", this_entry_year))
    next
  }
  hosps[,visit_num := as.character(visit_num)]
  hosps <- merge(hosps, codes, by = "visit_num")
  
  
  sums <- hosps[,lapply(.SD, sum), by = .(QID, AYEAR), 
                .SDcols = names(codes)[2:ncol(codes)]]
  
  write_fst(sums, paste0("./data/cache_data/hosps_by_qid/medicare_hosp_counts_entry_", this_entry_year, ".fst"))
    
}

end <- Sys.time(); end - start

