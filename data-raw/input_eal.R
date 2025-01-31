## code to prepare `input_eal` dataset goes here
ct <- cols(
  system = col_character(),
  num_stories = col_double(),
  design_s = col_character(),
  design_ns = col_character(),
  model = col_character(),
  intervention = col_double(),
  loss_ratio = col_double(),
  repair_cost = col_double(),
  re_occupancy_time = col_double(),
  functional_recovery_time = col_double()
)

input_eal <- readr::read_csv("data-raw/input_eal.csv", col_types = ct)

usethis::use_data(input_eal, overwrite = TRUE)


vars <- c(
  system = "Structural System (Lateral Force Resisting System)",
  num_stories = "Number of stories",
  design_s = "Structural design intervention (baseline/RC IV/backup-frame)",
  design_ns = "Nonstructural design intervention (baseline/nsfr)",
  model = "Model name (system-num_stories-design_s-design_ns)",
  intervention = "Recovery-based design intervention (0/1)",
  loss_ratio = "Loss ratio",
  repair_cost = "Annualized repair cost (dollars)",
  re_occupancy_time = "Annualized re-occupancy time (days)",
  functional_recovery_time = "Annualized functional recovery time (days)"
)

glue::glue("#'   \\item{[colname]}{[coldesc]}",
           colname = names(vars),
           coldesc = vars,
           .open = "[",
           .close = "]")
