data <- read.table("data-raw/data.txt")
usethis::use_data(data, overwrite = T)
