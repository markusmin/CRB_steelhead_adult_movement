# 02_clean_detection_histories: Run files 9-13 through script using hyak

source("R/functions/create_detection_history.R")

file_paths_3 <- c("Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_9_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_10_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_11_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_12_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_13_2025-01-28.csv")

output_3 <- create_detection_history(file_paths = file_paths_3)
write.csv(output_3$det_hist, "intermediate_outputs/CTH_9_13_complete_det_hist.csv", row.names = F)
write.csv(output_3$event_site_metadata, "intermediate_outputs/CTH_9_13_event_site_metadata.csv", row.names = F)
