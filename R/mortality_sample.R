#' Dataset giving exposure-to-death (population size) and number of deaths for males in eleven european countries,
#' at males with ages ranging from 30 years old to 85 years old
#'
#' @format Two different lists related to the reduced (subsample) population size and reduced number of deaths in
#'         eleven european countries, for male people aged 30 years old to 85 years old between 1908 and 2020.
#'         The data were exported from the Human Mortality Database (HMD).
#'
#' An evolving data frame of exposure-to-death and number of deaths in Belgium, Switzerland, Denmark, Spain, Finland,
#' France, United Kingdom, Italia, The Netherlands, Norway and Sweden.
#' \describe{
#'   \item{XP_red}{A list of nine elements (one for each country) giving a subset of the exposure-to-death (or reduced population size),
#'                 each element having 56 rows (ages 30-85) and 113 columns (period 1908-2020)}
#'   \item{DX_red}{A list of nine elements (one for each country) giving a subset of the number of deaths, each element having 56 rows
#'                 (ages 30-85) and 113 columns (period 1908-2020)}
#'   \item{country_names}{A list of nine elements giving the names of the countries, in the same order as the elements in other lists}
#' }
#' @source \url{www.mortality.org}
"mortality_sample"
