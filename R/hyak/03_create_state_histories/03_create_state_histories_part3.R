# 03_create_detection_histories_part3: Create the full state history for the third quartile of the data

source("R/functions/create_state_history.R")

output <- create_state_history(data_quartile = 3)
write.csv(output, "intermediate_outputs/states_complete_part3.csv", row.names = F)