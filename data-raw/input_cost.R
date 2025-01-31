## code to prepare `input_cost` dataset goes here
ct <- cols(
  system = col_character(),
  num_stories = col_double(),
  design_s = col_character(),
  design_ns = col_character(),
  model = col_character(),
  intervention = col_double(),
  c_s = col_double(),
  c_ns = col_double(),
  construction = col_double(),
  project = col_double()
)

input_cost <- readr::read_csv("data-raw/input_cost.csv", col_types=ct)

usethis::use_data(input_cost, overwrite = TRUE)

vars <- c(
  system = "Structural system (Lateral Force Resisting System)",
  num_stories = "Number of stories",
  design_s = "Structural design intervention (baseline/RC IV/backup-frame)",
  design_ns = "Nonstructural design intervention (baseline/nsfr)",
  model = "Model name (system-num_stories-design_s-design_ns)",
  intervention = "Recovery-based design intervention (0/1)",
  c_s = "Structural construction cost",
  c_ns = "Nonstructural construction cost",
  construction = "Total construction cost",
  project = "Total project cost"
)

glue::glue("#'   \\item{[colname]}{[coldesc]}",
           colname = names(vars),
           coldesc = vars,
           .open = "[",
           .close = "]")
