#' Deaths statistics in 11 european countries
#'
#' @format Dataset providing the exposure-to-death (population size) and number of deaths for males in 11 european countries,
#'         between 1908 and 2020, with ages ranging from 30 years old to 85 years old. Exported from the Human Mortality Database (HMD).
#'         The two first lists relate to some subsample of the population size and number of deaths in those countries, with
#'         random sampling from the original dataset.
#'
#' An evolving data frame of exposure-to-death and number of deaths in Belgium, Switzerland, Denmark, Spain, Finland,
#' France, United Kingdom, Italia, The Netherlands, Norway and Sweden.
#' \describe{
#'   \item{XP}{A list of eleven elements (one for each country) giving a subset of the exposure-to-death (or reduced population size),
#'                 each element having 56 rows (ages 30-85) and 113 columns (period 1908-2020)}
#'   \item{DX}{A list of eleven elements (one for each country) giving a subset of the number of deaths, each element having 56 rows
#'                 (ages 30-85) and 113 columns (period 1908-2020)}
#'   \item{names}{A list of eleven elements giving the names of the countries, in the same order as the elements in other lists}
#' }
#' @source \url{https://www.mortality.org}
"mortality_sample"
