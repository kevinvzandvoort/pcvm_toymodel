uganda_data = fread("./data/carriage_uganda.csv")
uganda_data[, c("S", "VT", "NVT") := .(N-TOTALST, VT10, TOTALST-VT10)]
uganda_data = setAgeBreaks(c(0, 2, 5, 15)) %>% cbind(uganda_data[, -"age"])
data_unvacll = uganda_data = uganda_data[, -c("VT10", "VT13", "TOTALST", "N")]

#source("./data/data_case_carrier_ratios.R")
source("./data/data_contact_matrix_ethiopia.R")
