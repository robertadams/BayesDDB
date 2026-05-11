# Setup dependencies
usethis::use_package("dplyr", type = "Imports")
usethis::use_package("overlapping", type = "Imports")
usethis::use_package("ggplot2", type = "Imports")
usethis::use_package("ggpubr", type = "Imports")

usethis::use_package("knitr", type = "Suggests")
usethis::use_package("rmarkdown", type = "Suggests")
usethis::use_package("rjags", type = "Suggests")
usethis::use_package("testthat", min_version = "3.0.0", type = "Suggests")

# Check URLS
urlchecker::url_check() # Currently fails due to private GitHub repo addresses

# Update wordlist
usethis::use_spell_check()
spelling::update_wordlist()

# Add GitHub actions checks
usethis::use_github_action()

# Jarl GHA
# usethis::use_github_action(url = "https://github.com/etiennebacher/setup-jarl/blob/main/examples/jarl-check.yml")

# Air GHA
# usethis::use_github_action(url = "https://github.com/posit-dev/setup-air/blob/main/examples/format-check.yaml")

# Checks
devtools::load_all()
devtools::document()
devtools::test()
devtools::check()
devtools::check(remote = TRUE, manual = TRUE)

# Check against CRAN's win builder service
# By default will be sent to maintainer, but can be overriden setting email arg and commenting out maintainer in DESCRIPTION file
# devtools::check_win_devel(email = "")
# devtools::check_win_release(email = "")
