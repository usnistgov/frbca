
<!-- README.md is generated from README.Rmd. Please edit that file -->

# frbca

<!-- badges: start -->
<!-- badges: end -->

The goal of frbca is to conduct benefit-cost analysis (BCA) for
buildings archetypes designed for functional recovery. The inputs
include:

- Expected Annualized Losses (EALs) for (i) a baseline (code-conforming)
  design and (ii) recovery-based design interventions
- Structural and nonstructural costs for both (i) and (ii)
- A list of analysis parameters, including discount rate, time horizon,
  and economic analysis parameters such as business income

## Installation

You can install the development version of frbca like so:

You can install frbca from [GitHub](https://github.com) with:

``` r
devtools::install_github("juanfung/frbca")
```

## Expected inputs

To run the analysis, provide two inputs in tabular (eg, `data.frame`)
form:

### Expected annualized losses (EALs) from a performance assessment

Each row is an archetype model. Expected (minimum) columns:

    | variable                 | type  | description                        |
    |--------------------------|-------|------------------------------------|
    | model                    | <chr> | unique model name or id            |
    | intervention             | {0,1} | baseline design = 0; otherwise = 1 |
    | num_stories              | <dbl> | number of stories                  |
    | loss_ratio               | <dbl> | loss ratio                         |
    | repair_costs             | <dbl> | repair costs, in dollars           |
    | re_occupancy_time        | <dbl> | re-occupancy time, in days         |
    | functional_recovery_time | <dbl> | functional recovery time, in days  |

### Construction costs

Each row is an archetype model. Expected (minimum) columns:

    | variable     | type  | description                                  |
    |--------------|-------|----------------------------------------------|
    | model        | <chr> | unique model name or id                      |
    | intervention | {0,1} | baseline design = 0; otherwise = 1           |
    | num_stories  | <dbl> | number of stories                            |
    | c_s          | <dbl> | structural construction costs, in dollars    |
    | c_ns         | <dbl> | nonstructural construction costs, in dollars |

In addition, several parameters are required in list form. The following
are the base parameters required under sub-list
`input_param$parameters$base`:

    | parameter    | type    | description                                        |
    |--------------|---------|----------------------------------------------------|
    | floor_area   | <dbl>   | total square footage per story                     |
    | delta        | <delta> | discount rate                                      |
    | T            | <dbl>   | planning horizon, in years                         |
    | bi           | <dbl>   | business income                                    |
    | ri           | <dbl>   | rental income                                      |
    | displacement | <dbl>   | occupant-incurred costs of displacement            |
    | tenant       | <dbl>   | number of tenants per square foot                  |
    | recapture    | <dbl>   | recapture rate for rental income                   |
    | sc           | <sc>    | supply-chain multiplier for business income losses |

In addition, parameter values for sensitivity analysis should be
provided under the sub-list `input_param$parameters$sensitivity`. For
example, sensitivity analysis for `delta` would be conducted based on
provided values `input_param$parameters$sensitivity$delta$low` and
`input_param$parameters$sensitivity$delta$high`.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
## NOT RUN
library(frbca)

## assuming input data is loaded and parameters are set...

## run analysis
output <- frbca::frbca(input_eal, input_cost, input_param)

## generate results plot with sensitivity analysis
frbca::plot_frbca(output, n_floors=4, system='RCMF')
```
