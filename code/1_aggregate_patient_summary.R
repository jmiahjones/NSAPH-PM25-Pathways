## Aggregate Medicare Denominator By Zip

library(NSAPHutils)
library(NSAPHplatform)
set_threads()
library(fst)
library(data.table)

for (year in 1999:2016) {
  x <- read_data("../data/denom/", years = year, columns = c('zip',
                                                             'year',
                                                             'sex',
                                                             'race',
                                                             'age',
                                                             'statecode',
                                                             'dead',
                                                             'pm25_ensemble',
                                                             'poverty',
                                                             'popdensity',
                                                             'medianhousevalue',
                                                             'pct_blk',
                                                             'medhouseholdincome',
                                                             'pct_owner_occ',
                                                             'hispanic',
                                                             'education',
                                                             'population'))
  x[, m_count := sum(sex == 1), by = "zip"]
  x[, f_count := sum(sex == 2), by = "zip"]
  x[, white_count := sum(race == 1), by = "zip"]
  x[, black_count := sum(race == 2), by = "zip"]
  x[, hispanic_count := sum(race == 5), by = "zip"]
  x[, asian_count := sum(race == 4), by = "zip"]
  x[, native_count := sum(race == 6), by = "zip"]
  x <- x[, .(statecode = max(statecode),
             year = max(year),
             total_count = .N,
             m_count = max(m_count),
             f_count = max(f_count),
             mean_age = mean(age, na.rm = T),
             white_count = max(white_count),
             black_count = max(black_count),
             hispanic_count = max(hispanic_count),
             asian_count = max(asian_count),
             native_count = max(native_count),
             deaths = sum(dead),
             pm25 = max(pm25_ensemble),
             poverty = max(poverty),
             popdensity = max(popdensity),
             medianhousevalue = max(medianhousevalue),
             pct_blk = max(pct_blk),
             medhouseholdincome = max(medhouseholdincome),
             pct_owner_occ = max(pct_owner_occ),
             hispanic = max(hispanic),
             education = max(education),
             population = max(population)),
         by = "zip"]
  write_fst(x, paste0("../data/by_zip/data/patient_summary_by_zip_",year, ".fst"))
}