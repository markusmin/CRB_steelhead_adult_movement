# 03_create_detection_histories_part2: Create the full state history for the second quartile of the data

source("R/functions/create_state_history.R")

output <- create_state_history(data_quartile = 2)
write.csv(output, "intermediate_outputs/states_complete_part2.csv", row.names = F)