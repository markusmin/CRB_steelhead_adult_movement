# 02_clean_detection_histories: Run files 5-8 through script using hyak

source("R/functions/create_detection_history.R")

file_paths_2 <- c("Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_5_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_6_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_7_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_8_2025-01-28.csv")

output_2 <- create_detection_history(file_paths = file_paths_2)
write.csv(output_2$det_hist, "intermediate_outputs/CTH_5_8_complete_det_hist.csv", row.names = F)
write.csv(output_2$event_site_metadata, "intermediate_outputs/CTH_5_8_event_site_metadata.csv", row.names = F)
