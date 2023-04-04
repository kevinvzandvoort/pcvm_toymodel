pacman::p_load(magrittr, data.table)

#Estimate
HARGEISA_POPULATION = 1200000

#' Get contact data from Prem et al (Ethiopia)
contact_data = readRDS("./data/eth_contacts.RDS")
ethiopia_contact_matrix_agegroups = seq(0, 75, 5) %>% set_units("years") %>% setAgeBreaks(maxage = set_units(120, "years"))
colnames(contact_data) = ethiopia_contact_matrix_agegroups[, name]
contact_data = as.data.table(contact_data)
contact_data[, contactor_age_group := ethiopia_contact_matrix_agegroups[, name]]
ethiopia_contact_matrix = contact_data %>% melt(id.vars = "contactor_age_group", variable.name = "contactee_age_group")
ethiopia_contact_matrix[, c("contactor_age_group", "contactee_age_group") :=
                          .(factor(contactor_age_group, ethiopia_contact_matrix_agegroups[, name]),
                            factor(contactee_age_group, ethiopia_contact_matrix_agegroups[, name]))]

pacman::p_load(magrittr, data.table)

#' Get data from UNWPP for Somalia
population_data = fread("./data/somalia_popdata.csv")
population_data = population_data %>% melt(id.vars = c("country", "year"), variable.name = "age")
population_data[, value := value * 1000]
population_data[, age := as.numeric(as.character(age))]
population_data = population_data[country == "Somalia" & year == 2020] %>%
  cbind(c(c(1:100) %>% set_units("years")) %>% setAgeBreaks(maxage = set_units(120, "years")), .)
hargeisa_population_data = population_data[, -c("country", "year", "age")]
hargeisa_population_data[, value := value * HARGEISA_POPULATION/sum(value)]
hargeisa_contact_matrix_agegroups = ethiopia_contact_matrix_agegroups
