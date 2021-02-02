#link hosp, denom counts

library(NSAPHutils)
library(NSAPHplatform)
set_threads()
library(fst)
library(data.table)

for (year in 2000:2016) {
  denom <- read_data("../data/by_zip/data/", years = year, 
                     columns = fst.metadata("../data/by_zip/data/patient_summary_by_zip_2000.fst")$columnNames)
  denom <- denom[statecode %in% c("NY","MA","CT","RI","NH","ME","VT")]
  hosps <- read_data("../data/cache_data/", years = year,
                     columns = fst.metadata("../data/cache_data/medicare_hosp_counts_2000.fst")$columnNames)
  out <- merge(denom, hosps, by = "zip", all.x = T)
  out[, zip := int_to_zip_str(zip)]
  for (var in names(hosps)[2:285]) {
    out[is.na(get(var)), (var) := 0]
  }
  
  write_fst(out, paste0("../data/out_dir/ne_ccs_hosp_counts_", year, ".fst"))
}