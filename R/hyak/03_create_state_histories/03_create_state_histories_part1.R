# 03_create_detection_histories_part1: Create the full state history for the first quartile of the data

source("R/functions/create_state_history.R")

output <- create_state_history(data_quartile = 1)
write.csv(output, "intermediate_outputs/states_complete_part1.csv", row.names = F)