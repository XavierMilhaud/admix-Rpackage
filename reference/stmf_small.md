# Dataset of Short-term Mortality Fluctuations (STMF) from HMD

Restricted to 6 countries: Belgium, France, Italy, Netherlands, Spain,
Germany. Weekly death counts provide the most objective and comparable
way of assessing the scale of short-term mortality elevations across
countries and time. Extraction date from the Human Mortality Database
(HMD): 09/21/2020.

## Usage

``` r
stmf_small
```

## Format

A data frame with 88146 rows and 19 variables:

- CountryCode:

  Mortality database country code

- Year:

  Year

- Week:

  Week number

- Sex:

  Gender ('m': male, 'f': female, 'b': both)

- D0_14:

  Age range 0-14

- D15_64:

  Age range 15-64

- D65_74:

  Age range 65-74

- D75_84:

  Age range 75-84

- D85p:

  Age range 85-+

- DTotal:

  Count of deaths for all ages combined

- R0_14:

  Crude death rate for age range 0-14

- R15_64:

  Crude death rate for age range 15-64

- R65_74:

  Crude death rate for age range 65-74

- R75_84:

  Crude death rate for age range 75-84

- R85p:

  Crude death rate for age range 85-+

- RTotal:

  Crude death rate for all ages combined

- Split:

  Indicates if data were split from aggregated age groups (0 if the
  original data has necessary detailed age scale). For example, if the
  original age scale was 0-4, 5-29, 30-65, 65+, then split will be equal
  to 1

- SplitSex:

  Indicates if the original data are available by sex (0) or data are
  interpolated (1)

- Forecast:

  Equals 1 for all years where forecasted population exposures were used
  to calculate weekly death rates

## Source

<https://www.mortality.org>
