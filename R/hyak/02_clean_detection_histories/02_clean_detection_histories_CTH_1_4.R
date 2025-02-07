# 02_clean_detection_histories: Run files 1-4 through script using hyak

source("R/functions/create_detection_history.R")

file_paths_1 <- c("Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_1_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_2_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_3_2025-01-28.csv",
                  "Data/PTAGIS_queries/complete_tag_histories/CTH_tag_codes_4_2025-01-28.csv")

output_1 <- create_detection_history(file_paths = file_paths_1)
write.csv(output_1$det_hist, "intermediate_outputs/CTH_1_4_complete_det_hist.csv", row.names = F)
write.csv(output_1$event_site_metadata, "intermediate_outputs/CTH_1_4_event_site_metadata.csv", row.names = F)