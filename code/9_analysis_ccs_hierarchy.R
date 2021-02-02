library(dplyr)
ccs_hierarchy_tablea <- read.csv("./data/ccs-user-guide-table-a.csv", header=T, 
                                 sep=",", quote='"', stringsAsFactors = F,
                                 col.names = c("ccs", "ICD.Chapter"))

translation_df <- lapply(seq_along(ccs_hierarchy_tablea$CCS.Codes), function(idx){
  x <- ccs_hierarchy_tablea$CCS.Codes[idx]
  name <- ccs_hierarchy_tablea$ICD.Chapter[idx]
  data.frame(name=name, ccs=eval(parse(text=paste0("c(", gsub("-", ":", x), ")")))
  )
}) %>% dplyr::bind_rows() %>%
  mutate(ccs = paste0("ccs_", ccs))
