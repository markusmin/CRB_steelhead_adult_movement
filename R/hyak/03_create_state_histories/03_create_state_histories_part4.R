# 03_create_detection_histories_part4: Create the full state history for the fourth quartile of the data

source("R/functions/create_state_history.R")

output <- create_state_history(data_quartile = 4)
write.csv(output, "intermediate_outputs/states_complete_part4.csv", row.names = F)